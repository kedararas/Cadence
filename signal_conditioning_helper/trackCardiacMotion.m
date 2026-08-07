function [videoTracked, D, videoNorm] = trackCardiacMotion(video, varargin)
%TRACKCARDIACMOTION  Marker-free non-rigid motion tracking for optical mapping.
%
%   [videoTracked, D, videoNorm] = trackCardiacMotion(VIDEO)
%
%   Implements the two-stage marker-free motion-tracking scheme used in
%   Christoph & Luther (2018) and Kappadan et al. (J Physiol 2023):
%
%     (1) Locally normalize / contrast-enhance every frame so the tracker
%         keys on tissue texture (trabeculae, vessels, dye speckle) rather
%         than the low-spatial-frequency action-potential wavefront.
%         Local mean is subtracted and the residual divided by a local
%         standard deviation; this is a high-pass + local-variance
%         equalization that flattens the slow fluorescence changes caused
%         by the AP while preserving the structural texture that actually
%         moves with the tissue.
%
%     (2) Compute a dense, non-rigid displacement field from each test
%         frame to a single reference frame using multi-level diffeomorphic
%         demons registration (IMREGDEMONS). The field is computed on the
%         NORMALIZED video but applied (with IMWARP) to the ORIGINAL
%         fluorescence frames, so downstream electrophysiology analysis
%         sees motion-corrected raw data with its fluorescence content
%         untouched.
%
%   Inputs
%     VIDEO    [H x W x T]  fluorescence stack (uint16, single, or double)
%
%   Name / Value options
%     'RefFrame'       frame index (or [] for auto)            default 1
%                      Use the frame just before a pacing stimulus for
%                      paced data; any frame is fine for VF.
%     'SigmaLP'        Gaussian sigma (px) for local mean      default 8
%     'SigmaStd'       Gaussian sigma (px) for local std       default 8
%     'ClipStd'        clip normalized frame to +/- this value default 3
%     'PyramidLevels'  demons multi-resolution levels          default 3
%     'Iterations'     demons iterations per level             default [100 50 25]
%     'Smoothing'      demons elastic regularization (px)      default 1.5
%     'Mask'           [H x W] logical tissue mask             default []
%                      Outside-mask pixels are set to 0 in the normalized
%                      video so they contribute nothing to the registration.
%     'ReturnNorm'     also return normalized stack (large)    default false
%     'ReturnDisp'     also return displacement field D (large) default false
%                      D is [H x W x 2 x T] single — ~262 MB for 128x128x2000.
%                      Leave false (default) unless you need the raw field.
%     'Verbose'        progress output                         default true
%
%   Outputs
%     videoTracked  [H x W x T]  motion-corrected original fluorescence
%     D             [H x W x 2 x T]  per-frame displacement field (px), or
%                   [] when 'ReturnDisp' is false (default).
%                   D(:,:,1,t) is the x (column) displacement,
%                   D(:,:,2,t) is the y (row) displacement, such that
%                   IMWARP(testFrame, D(:,:,:,t)) approximates the
%                   reference frame.
%     videoNorm     [H x W x T]  (only if 'ReturnNorm', true) the locally
%                   normalized stack the tracker actually operated on.
%
%   Example
%     load('contracting_heart.mat');    % v: [H x W x T] OAP stack, fs
%     mask = v(:,:,1) > graythresh(v(:,:,1))*double(max(v(:)));
%     [vT, D] = trackCardiacMotion(v, 'RefFrame', 1, 'Mask', mask);
%     % Ensemble average at a pixel:
%     trace = squeeze(vT(80, 64, :));
%     [ap, M, t] = ensembleAvgAP(trace, fs, 'Plot', true);
%
%   Notes
%     * Requires Image Processing Toolbox (imgaussfilt, imregdemons, imwarp).
%     * Memory: D is 4-D and can be very large; for long recordings consider
%       chunking the video or discarding D after warping.
%     * For out-of-plane motion, 2-D tracking is fundamentally limited; this
%       routine, like the published method, corrects in-plane motion only.
%
% ------------------------------------------------------------------------

    % ---------- parse inputs ----------
    p = inputParser;
    p.addRequired('video', @(x) ndims(x)==3 && isnumeric(x));
    p.addParameter('RefFrame',       1);
    p.addParameter('SigmaLP',        8,    @(x) x>0);
    p.addParameter('SigmaStd',       8,    @(x) x>0);
    p.addParameter('ClipStd',        3,    @(x) x>0);
    p.addParameter('PyramidLevels',  3,    @(x) x>=1 && x==round(x));
    p.addParameter('Iterations',     [100 50 25]);
    p.addParameter('Smoothing',      1.5,  @(x) x>=0);
    p.addParameter('Mask',           []);
    p.addParameter('ReturnNorm',     false, @islogical);
    p.addParameter('ReturnDisp',     false, @islogical);
    p.addParameter('Verbose',        true,  @islogical);
    p.parse(video, varargin{:});
    o = p.Results;

    video = single(video);
    [H, W, T] = size(video);

    % ---------- blank-frame guard ----------
    % Acquisition stacks often contain blank (all-zero / near-constant) frames
    % — e.g. a leading warm-up frame or dropped frames. Registering against a
    % blank reference, or registering a blank frame, yields a garbage
    % displacement field. Flag them up front (a frame is "blank" if its spatial
    % std is a tiny fraction of the recording's typical frame std).
    frameStd  = squeeze(std(reshape(video, H*W, T), 0, 1));   % [T x 1]
    typicalStd = median(frameStd(frameStd > 0));
    if isempty(typicalStd) || ~isfinite(typicalStd)
        error('trackCardiacMotion:allBlank', 'All frames are blank/constant.');
    end
    isBlank = frameStd < 1e-3 * typicalStd;
    if o.Verbose && any(isBlank)
        fprintf('trackCardiacMotion: %d blank frame(s) detected (e.g. %s) — passed through unchanged.\n', ...
                nnz(isBlank), mat2str(find(isBlank, 5)'));
    end

    if isempty(o.RefFrame)
        % Auto: pick the frame with minimum total gradient change relative
        % to its neighbours (typically a diastolic frame).
        o.RefFrame = autoPickReference(video, isBlank);
    end
    assert(o.RefFrame >= 1 && o.RefFrame <= T, ...
           'RefFrame out of bounds.');

    % If the requested reference is blank, fall back to an auto-picked one.
    if isBlank(o.RefFrame)
        newRef = autoPickReference(video, isBlank);
        warning('trackCardiacMotion:blankRef', ...
            'RefFrame %d is blank; using auto-picked reference frame %d instead.', ...
            o.RefFrame, newRef);
        o.RefFrame = newRef;
    end

    % ---------- tissue mask ----------
    if isempty(o.Mask)
        mask = true(H, W);
    else
        mask = logical(o.Mask);
        assert(isequal(size(mask), [H, W]), 'Mask must be H x W.');
    end

    % ---------- build normalized reference ----------
    refNorm = localNormalize(video(:,:,o.RefFrame), ...
                             o.SigmaLP, o.SigmaStd, o.ClipStd);
    refNorm(~mask) = 0;

    % ---------- pre-allocate outputs ----------
    videoTracked = zeros(H, W, T, 'single');
    % D is only allocated when the caller explicitly requests it.
    % For a 128x128x2000 recording this saves ~262 MB.
    if o.ReturnDisp
        D = zeros(H, W, 2, T, 'single');
    else
        D = [];
    end
    if o.ReturnNorm
        videoNorm = zeros(H, W, T, 'single');
    else
        videoNorm = [];
    end

    % ---------- main loop ----------
    if o.Verbose
        fprintf('trackCardiacMotion: %d frames, ref=%d, %dx%d\n', ...
                T, o.RefFrame, H, W);
        tic;
    end

    reportStep = max(1, round(T / 20));   % pre-compute, not per iteration

    for k = 1:T
        % Blank frame: nothing to register against. Pass it through unchanged
        % with an identity field rather than registering noise-to-reference.
        if isBlank(k)
            videoTracked(:,:,k) = video(:,:,k);
            % D(:,:,:,k) and videoNorm(:,:,k) stay at their zero-initialised
            % values, which is the correct identity / blank result.
            continue;
        end

        % (1) local normalization on THIS frame
        testNorm = localNormalize(video(:,:,k), ...
                                  o.SigmaLP, o.SigmaStd, o.ClipStd);
        testNorm(~mask) = 0;

        if o.ReturnNorm
            videoNorm(:,:,k) = testNorm;
        end

        % (2) non-rigid registration: testNorm --> refNorm
        if k == o.RefFrame
            Dk = zeros(H, W, 2, 'single');       % identity displacement
        else
            Dk = imregdemons(testNorm, refNorm, o.Iterations, ...
                'PyramidLevels',        o.PyramidLevels, ...
                'AccumulatedFieldSmoothing', o.Smoothing, ...
                'DisplayWaitbar',       false);
        end
        if o.ReturnDisp
            D(:,:,:,k) = Dk;
        end

        % (3) warp the ORIGINAL (un-normalized) frame with Dk
        videoTracked(:,:,k) = imwarp(video(:,:,k), Dk, ...
                                     'FillValues', 0, 'SmoothEdges', true);

        if o.Verbose && mod(k, reportStep) == 0
            fprintf('  %4d / %4d  (%.1f%%)  elapsed %.1fs\n', ...
                    k, T, 100*k/T, toc);
        end
    end

    if o.Verbose
        fprintf('Done in %.1f s.\n', toc);
    end
end

% ========================================================================
function In = localNormalize(I, sigmaLP, sigmaStd, clipStd)
%LOCALNORMALIZE  High-pass + local-variance equalization.
%
%   In = (I - mu_local) / (sigma_local + eps), clipped to +/- clipStd.
%
%   The Gaussian low-pass removes the slowly-varying component — the AP
%   itself is a ~cm-scale wave of uniform fluorescence drop, so it sits
%   almost entirely in the low-pass band. The residual is divided by a
%   local standard deviation so that dim and bright regions of the heart
%   contribute comparable signal to the optical-flow cost function.
%   The tracker therefore locks onto texture (trabeculae, vessels, dye
%   speckle) that co-moves with the tissue, not onto the electrical wave.

    I = single(I);

    mu    = imgaussfilt(I, sigmaLP);
    dI    = I - mu;
    sigma = sqrt(imgaussfilt(dI.^2, sigmaStd) + eps('single'));
    In    = dI ./ sigma;

    In(In >  clipStd) =  clipStd;
    In(In < -clipStd) = -clipStd;
end

% ========================================================================
function idx = autoPickReference(video, isBlank)
%AUTOPICKREFERENCE  Heuristic: choose the frame whose spatial gradient
% magnitude is closest to the median across time. In a periodically paced
% recording this tends to land on a diastolic frame.
%
% Blank frames (isBlank) are excluded from both the median and the choice so
% the reference is never a zero/constant frame.
%
% Vectorised: forward differences along x and y are computed across the
% full 3D stack in two calls, avoiding a per-frame loop.
    video = single(video);
    T = size(video, 3);
    if nargin < 2 || isempty(isBlank)
        isBlank = false(T, 1);
    end
    gx = diff(video, 1, 2);          % [H x W-1 x T]  forward diff along cols
    gy = diff(video, 1, 1);          % [H-1 x W x T]  forward diff along rows
    gxm = squeeze(mean(mean(gx.^2, 1), 2));   % [T x 1]  mean sq-grad per frame
    gym = squeeze(mean(mean(gy.^2, 1), 2));   % [T x 1]
    g   = sqrt(gxm + gym);

    % Score only valid frames; force blanks to never win.
    valid = ~isBlank(:);
    medG  = median(g(valid));
    score = abs(g - medG);
    score(~valid) = inf;
    [~, idx] = min(score);
end
