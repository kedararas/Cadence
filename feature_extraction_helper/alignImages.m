function [movingReg, tform] = alignImages(moving, fixed, varargin)
%ALIGNIMAGES  Register a moving cardiac optical-mapping image to a fixed
%             reference using intensity-based or feature-based registration.
%
%   [movingReg, tform] = alignImages(moving, fixed)
%   [movingReg, tform] = alignImages(moving, fixed, Name, Value, ...)
%
%   Designed for dual-camera cardiac optical mapping: both cameras image
%   the same dye-stained tissue (monomodal). Both cameras share the same
%   pixel size, so a RIGID transform (translation + rotation, NO scale)
%   is the correct default. 'similarity' is available but can produce a
%   degenerate squeeze where the optimizer shrinks tissue into the dark
%   background — only use it when the cameras have different magnifications.
%
%   Inputs
%     moving, fixed   Grayscale or RGB images, any numeric class (uint8,
%                     uint16, single, double).
%
%   Name-Value options
%     'Method'      'intensity' (default) | 'feature'
%                   Intensity-based is recommended for cardiac OM images
%                   (soft texture; few reliable keypoints for SURF).
%     'Transform'   'rigid' (default) | 'similarity' | 'affine'
%                   Rigid = translation + rotation only. Correct for
%                   dual cameras with the same optics. Use 'similarity'
%                   only when cameras have genuinely different zoom levels
%                   (it adds a uniform scale factor that can diverge).
%     'ShowPlot'    true (default) | false
%
%   Outputs
%     movingReg   Moving image warped into the fixed image coordinate frame.
%     tform       Geometric transform object (affine2d or affinetform2d
%                 depending on MATLAB version). Pass to imwarp/
%                 transformPointsForward to align full video stacks.
%
%   Example — align two camera images and apply to full video:
%     [cam2reg, tform] = alignImages(cam2_image, cam1_image);
%     outRef = imref2d(size(cam1_image));
%     for k = 1:T
%         cam2_video_reg(:,:,k) = imwarp(cam2_video(:,:,k), tform, ...
%                                        'OutputView', outRef, 'FillValues', 0);
%     end

% ── Parse inputs ─────────────────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'Method',    'intensity', ...
    @(s) any(strcmpi(s, {'intensity','feature'})));
addParameter(p, 'Transform', 'rigid', ...
    @(s) any(strcmpi(s, {'rigid','similarity','affine'})));
addParameter(p, 'ShowPlot',  true, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});
opt           = p.Results;
opt.Method    = lower(opt.Method);
opt.Transform = lower(opt.Transform);
showPlot      = logical(opt.ShowPlot);

% ── Step 1: Convert to single-precision grayscale in [0, 1] ──────────────────
% The previous code did fixedG = fixed; movingG = moving;
% without any conversion. detectSURFFeatures and imregtform both require
% grayscale input; uint16 is not accepted by either function.
fixedG  = toGraySingle(fixed);
movingG = toGraySingle(moving);

% ── Step 2: Tissue mask — exclude dark background ────────────────────────────
% Cardiac OM images are a bright heart on a near-black background. The large
% uniform background region dominates the cost function and misleads the
% optimizer toward the trivial (zero-motion) solution.
% We zero out the background so registration locks onto the tissue only.
fixedMask  = makeTissueMask(fixedG);
movingMask = makeTissueMask(movingG);

% ── Step 3: Enhance local contrast (CLAHE) ───────────────────────────────────
% Cardiac texture (trabeculae, vessels, dye speckle) is low-contrast.
% CLAHE amplifies local texture within the tissue region, giving the
% optimizer (or feature detector) stronger gradient signal to lock onto.
fixedE  = enhanceTissue(fixedG,  fixedMask);
movingE = enhanceTissue(movingG, movingMask);

% ── Step 4: Register ─────────────────────────────────────────────────────────
switch opt.Method

    % ── Intensity-based (recommended for cardiac OM) ──────────────────────
    case 'intensity'
        [optimizer, metric] = imregconfig('monomodal');

        % Tuned for cardiac OM images:
        %  - More iterations: soft tissue gradients are shallow
        %  - Step 0.05: small enough to not overshoot the cost minimum
        %  - Relaxation: stabilises gradient-descent trajectory
        optimizer.MaximumIterations = 400;
        optimizer.MaximumStepLength = 0.1;
        optimizer.MinimumStepLength = 1e-6;
        optimizer.RelaxationFactor  = 0.6;

        % PyramidLevels=3 is safer than 5:
        %   5 levels on a 200×200 image → coarsest is ~6×6 (no texture left)
        %   3 levels → coarsest is ~25×25 (enough gradient signal)
        tform = imregtform(movingE, fixedE, opt.Transform, optimizer, metric, ...
            'PyramidLevels', 3);

        % ── Sanity check for similarity: guard against scale divergence ───────
        % A correctly registered dual-camera pair should have scale ≈ 1.
        % If the optimizer produced scale < 0.7 or > 1.4 the tissue was
        % squeezed/expanded into a degenerate solution; fall back to rigid.
        if strcmpi(opt.Transform, 'similarity')
            try
                if isprop(tform, 'A')          % R2022b+ affinetform2d
                    sc = sqrt(tform.A(1,1)^2 + tform.A(2,1)^2);
                else                           % pre-R2022b affine2d (.T matrix)
                    sc = sqrt(tform.T(1,1)^2 + tform.T(1,2)^2);
                end
                if sc < 0.7 || sc > 1.4
                    warning('alignImages:badScale', ...
                        ['Similarity fit produced scale=%.2f (expected ~1.0). ' ...
                         'The optimizer found a degenerate solution. ' ...
                         'Retrying with ''rigid'' to prevent squeeze.'], sc);
                    tform = imregtform(movingE, fixedE, 'rigid', optimizer, metric, ...
                        'PyramidLevels', 3);
                end
            catch
                % Could not inspect tform scale — continue with current result.
            end
        end

    % ── Feature-based (fallback / alternative) ────────────────────────────
    case 'feature'
        if ~exist('detectSURFFeatures', 'file')
            error('alignImages:noToolbox', ...
                ['Feature method requires Computer Vision Toolbox. ' ...
                 'Use ''Method'',''intensity'' instead.']);
        end

        % CRITICAL FIX: MetricThreshold was 500 — far too strict for the
        % soft texture of cardiac OM images. Lower to 50 and increase
        % NumOctaves so features are detected at multiple scales.
        ptsF = detectSURFFeatures(fixedE,  'MetricThreshold', 50, 'NumOctaves', 6);
        ptsM = detectSURFFeatures(movingE, 'MetricThreshold', 50, 'NumOctaves', 6);

        % Keep only the strongest points to reduce false matches
        ptsF = ptsF.selectStrongest(300);
        ptsM = ptsM.selectStrongest(300);

        [fF, vF] = extractFeatures(fixedE,  ptsF,  'Upright', false);
        [fM, vM] = extractFeatures(movingE, ptsM,  'Upright', false);

        % MaxRatio=0.75: looser than Lowe's 0.6 — cardiac images have
        % fewer distinctive keypoints so a stricter ratio drops too many.
        pairs = matchFeatures(fM, fF, 'Unique', true, 'MaxRatio', 0.75);

        if size(pairs, 1) < 4
            warning('alignImages:fewMatches', ...
                ['Only %d feature matches found — falling back to ' ...
                 'intensity-based registration.'], size(pairs, 1));
            [movingReg, tform] = alignImages(moving, fixed, ...
                'Method', 'intensity', ...
                'Transform', opt.Transform, ...
                'ShowPlot',  showPlot);
            return;
        end

        mM = vM(pairs(:,1));
        mF = vF(pairs(:,2));

        % estimateGeometricTransform2D was renamed in R2022b; handle both.
        try
            [tform, ~, ~] = estimateGeometricTransform2D( ...
                mM, mF, opt.Transform, 'MaxDistance', 6, 'Confidence', 99.5);
        catch
            [tform, ~, ~] = estgeotform2d( ...
                mM, mF, opt.Transform, 'MaxDistance', 6, 'Confidence', 99.5);
        end

        if showPlot
            figure('Name', 'Feature Matches', 'NumberTitle', 'off');
            showMatchedFeatures(fixedE, movingE, mF, mM, 'montage');
            title(sprintf('Feature matches used for %s registration', opt.Transform));
        end
end

% ── Step 5: Apply transform to original image ─────────────────────────────────
% The transform was computed on grayscale/enhanced images but is applied to
% the original (which may be uint16 or RGB) so downstream analysis sees the
% unmodified fluorescence intensities.
outRef    = imref2d(size(fixedG));
movingReg = imwarp(moving, tform, 'OutputView', outRef, 'FillValues', 0);

% ── Step 6: Visualise ─────────────────────────────────────────────────────────
if showPlot
    figure('Name', 'Alignment Result', 'NumberTitle', 'off', ...
           'Position', [100 100 1200 400]);

    subplot(1, 3, 1);
    imshow(fixedG, []);
    title('Fixed (CAM1)', 'FontWeight', 'bold');

    subplot(1, 3, 2);
    imshow(movingG, []);
    title('Moving (CAM2) — before', 'FontWeight', 'bold');

    subplot(1, 3, 3);
    imshowpair(fixedG, toGraySingle(movingReg), 'checkerboard');
    title(sprintf('Registered (%s, %s) — after', opt.Method, opt.Transform), ...
          'FontWeight', 'bold');

    sgtitle('Dual-Camera Cardiac OM Alignment');
    drawnow;
end
end


% ════════════════════════════════════════════════════════════════════════════
%  Local helpers
% ════════════════════════════════════════════════════════════════════════════

function I = toGraySingle(I)
%TOGRAYSINGLE  Convert any image class/colour-space to single [0,1] grayscale.
    % Colour → grayscale
    if size(I, 3) == 3
        I = rgb2gray(I);
    end
    % Any integer type → single
    I = single(I);
    % Normalise to [0, 1]
    lo = min(I(:));
    hi = max(I(:));
    if hi > lo
        I = (I - lo) / (hi - lo);
    else
        I = zeros(size(I), 'single');
    end
end


function mask = makeTissueMask(I)
%MAKETISSUEMASK  Segment heart tissue from background in a cardiac OM image.
%
%   Strategy: Otsu threshold (slightly relaxed to capture the full tissue
%   boundary), then morphological fill + close to handle the dark vessel
%   shadows inside the heart that would otherwise create holes.
    thresh = graythresh(I);
    % Relax to 60% of Otsu: cardiac images have a bright core but the edges
    % are dimmer than the automatic threshold expects.
    thresh = max(thresh * 0.6, 0.04);
    mask   = I > thresh;
    % Fill any interior holes (dark vessels, cannula shadows)
    mask   = imfill(mask, 'holes');
    % Close small boundary gaps
    mask   = imclose(mask, strel('disk', 8));
    % Dilate slightly to include the tissue boundary gradient
    mask   = imdilate(mask, strel('disk', 4));
end


function Ie = enhanceTissue(I, mask)
%ENHANCETISSUE  Apply CLAHE inside the tissue mask to boost local contrast.
%
%   CLAHE (Contrast-Limited Adaptive Histogram Equalisation) amplifies
%   the local texture of trabeculae and vessels so the registration
%   optimizer and SURF detector have sharper gradients to lock onto.
%   Background pixels are zeroed to prevent them from anchoring the fit.
    Ie = adapthisteq(I, 'ClipLimit', 0.03, 'NumTiles', [8 8]);
    Ie = Ie .* single(mask);    % zero background — tissue only
end
