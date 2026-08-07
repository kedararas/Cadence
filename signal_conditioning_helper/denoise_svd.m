function [denoised, info] = denoise_svd(data, K, varargin)
%DENOISE_SVD  Spatiotemporal SVD/PCA denoising of an optical-mapping stack.
%
%   denoised = denoise_svd(data)            % auto-select rank by energy
%   denoised = denoise_svd(data, K)         % keep K components
%   [denoised, info] = denoise_svd(data, K, Name, Value, ...)
%
%   The action-potential / Ca-transient signal is highly correlated across
%   space and beats, so it lives in a LOW-RANK subspace of the
%   [pixels x time] matrix.  Shot / read noise is spatially incoherent and
%   spreads across all the small singular values.  Keeping only the top
%   components reconstructs the coherent signal and discards most of the
%   noise WITHOUT the spatial blurring that spatial binning causes.
%
%   INPUTS
%     data   cmos data [rows x cols x frames]
%     K      number of singular components to keep.
%            []  or omitted -> auto (smallest K reaching 'EnergyThreshold')
%
%   Name/Value
%     'EnergyThreshold' fraction of variance to retain for auto-K (default 0.95)
%     'MaxRank'         cap on components considered/kept (default 15)
%     'MinComponents'   minimum K for the result to be considered usable
%                       (default 2).  If fewer components are selected, one
%                       mode dominates the variance and the reconstruction
%                       would be (near) rank-1 / spatially uniform, so SVD is
%                       NOT applied (see info.applied below).
%     'MaxPixelCorr'    max allowed mean inter-pixel correlation of the
%                       reconstruction (default 0.9995).  Even with K >=
%                       MinComponents the output can collapse to a near-uniform
%                       image; if the measured correlation reaches this value
%                       SVD is NOT applied.  Catches the spatial-collapse
%                       symptom directly, regardless of K.
%     'Mask'            [rows x cols] logical/numeric ROI; pixels outside are
%                       left unchanged and excluded from the decomposition
%
%   OUTPUTS
%     denoised  same size as data.  If SVD was not applied (info.applied ==
%               false), this is the INPUT data returned unchanged.
%     info      struct: .K (rank used), .applied (logical: true if SVD was
%               applied, false if the caller should fall back e.g. to spatial
%               binning), .singular_values, .energy_retained,
%               .energy_curve (cumulative energy vs rank)
%
%   METHOD
%     1. Reshape to M = [pixels x frames]; restrict to valid (in-mask,
%        finite) pixels.
%     2. Remove each pixel's temporal mean (denoise the dynamic part only).
%     3. Economy SVD via the top-MaxRank components (svds; falls back to a
%        full economy svd if svds is unavailable/fails).
%     4. Pick K (given or by cumulative singular-value energy).
%     5. Reconstruct with K components, add the per-pixel means back,
%        scatter into the original grid (out-of-mask pixels untouched).

    %% ---- parse inputs ----
    if nargin < 2, K = []; end
    p = inputParser;
    p.addParameter('EnergyThreshold', 0.95, @(v) isscalar(v) && v>0 && v<=1);
    p.addParameter('MaxRank',         15,   @(v) isscalar(v) && v>=1);
    p.addParameter('MinComponents',   2,    @(v) isscalar(v) && v>=1);
    p.addParameter('MaxPixelCorr',    0.9995,@(v) isscalar(v) && v>0 && v<=1);
    p.addParameter('Mask',            [],   @(v) isempty(v) || isnumeric(v) || islogical(v));
    p.parse(varargin{:});
    o = p.Results;

    data = double(data);
    [R, C, T] = size(data);

    %% ---- reshape to [pixels x frames] and select valid pixels ----
    M = reshape(data, R*C, T);                 % pixels x frames
    valid = all(isfinite(M), 2);
    if ~isempty(o.Mask)
        valid = valid & logical(o.Mask(:));
    end
    if ~any(valid)
        warning('denoise_svd:noValidPixels', 'No valid pixels; returning input.');
        denoised = data;  info = struct('K',0); return;
    end

    Mv = M(valid, :);                          % validPix x frames

    %% ---- remove per-pixel temporal mean ----
    mu = mean(Mv, 2);                          % validPix x 1
    Mc = Mv - mu;                              % centred

    %% ---- top-MaxRank SVD (economy) ----
    maxr = min([o.MaxRank, size(Mc,1)-1, size(Mc,2)-1]);
    maxr = max(maxr, 1);
    try
        [U, S, V] = svds(Mc, maxr);
        sv = diag(S);
    catch
        % Fallback: full economy SVD (slower, more memory) then truncate.
        [U, S, V] = svd(Mc, 'econ');
        sv = diag(S);
        keep = 1:min(maxr, numel(sv));
        U = U(:,keep); S = S(keep,keep); V = V(:,keep); sv = sv(keep);
    end

    %% ---- choose K ----
    % Normalise by the TRUE total variance (Frobenius norm of the centred
    % data), NOT by sum(sv.^2).  sv holds only the top-MaxRank singular
    % values, so sum(sv.^2) understates the total: a single dominant
    % component (e.g. photobleach drift, or a near-synchronous paced AP)
    % would then read as ~100% of the captured energy and force K=1 -- a
    % rank-1, spatially-uniform reconstruction (identical pixels everywhere).
    total_var = sum(Mc(:).^2);
    energy    = cumsum(sv.^2) / (total_var + eps);   % true cumulative fraction
    if isempty(K)
        Kuse = find(energy >= o.EnergyThreshold, 1, 'first');
        if isempty(Kuse), Kuse = numel(sv); end      % threshold not reached in MaxRank
    else
        Kuse = min(K, numel(sv));
    end
    %% ---- decide whether the decomposition is usable ----
    % Fewer than MinComponents means a single mode dominates the variance
    % (commonly photobleach drift on un-detrended data, or a near-synchronous
    % paced AP).  The reconstruction would then be (near) rank-1 / spatially
    % uniform, so flag applied=false and return the INPUT unchanged -- the
    % caller should fall back, e.g. to spatial binning + temporal filtering.
    applied = Kuse >= o.MinComponents;

    if ~applied
        warning('denoise_svd:notApplied', ...
            ['SVD not applied: only %d component(s) reached the energy ' ...
             'threshold (one mode dominates the variance). Returning the ' ...
             'input unchanged -- fall back to spatial binning + temporal ' ...
             'filtering. Detrend before SVD if you have not already.'], Kuse);
        denoised = data;                       % unchanged input
        info = struct( ...
            'K',               Kuse, ...
            'applied',         false, ...
            'singular_values', sv, ...
            'energy_retained', energy(min(max(Kuse,1), numel(energy))), ...
            'energy_curve',    energy, ...
            'pixel_corr',      NaN);
        fprintf('denoise_svd: NOT applied (K=%d < MinComponents=%d, %d valid pixels)\n', ...
                Kuse, o.MinComponents, sum(valid));
        return;
    end

    %% ---- reconstruct with K components ----
    Mc_hat = U(:,1:Kuse) * S(1:Kuse,1:Kuse) * V(:,1:Kuse)';
    Mv_hat = Mc_hat + mu;                      % add temporal means back

    %% ---- guard against a spatially-collapsed reconstruction ----
    % Even with K >= MinComponents the output can be near-uniform (every pixel
    % almost the same trace) if the retained modes are all spatially flat.
    % Measure the mean inter-pixel correlation on a sample of reconstructed
    % pixels; ~1 means the result collapsed and is not usable.
    pix_corr = sample_pixel_corr(Mv_hat);
    if pix_corr >= o.MaxPixelCorr
        applied = false;
        warning('denoise_svd:collapsed', ...
            ['SVD not applied: reconstruction is spatially collapsed (mean ' ...
             'inter-pixel correlation %.4f >= %.4f). Returning input unchanged ' ...
             '-- fall back to spatial binning + temporal filtering.'], ...
            pix_corr, o.MaxPixelCorr);
        denoised = data;                       % unchanged input
        info = struct( ...
            'K',               Kuse, ...
            'applied',         false, ...
            'singular_values', sv, ...
            'energy_retained', energy(min(Kuse, numel(energy))), ...
            'energy_curve',    energy, ...
            'pixel_corr',      pix_corr);
        fprintf('denoise_svd: NOT applied (collapsed, inter-pixel corr=%.4f)\n', pix_corr);
        return;
    end

    %% ---- scatter back to full grid (out-of-mask pixels untouched) ----
    Mout = M;
    Mout(valid, :) = Mv_hat;
    denoised = reshape(Mout, R, C, T);

    %% ---- info ----
    info = struct( ...
        'K',               Kuse, ...
        'applied',         true, ...
        'singular_values', sv, ...
        'energy_retained', energy(min(Kuse, numel(energy))), ...
        'energy_curve',    energy, ...
        'pixel_corr',      pix_corr);

    fprintf('denoise_svd: kept %d / %d components (%.1f%% var, %d px, inter-pixel corr=%.3f)\n', ...
            Kuse, numel(sv), 100*info.energy_retained, sum(valid), pix_corr);
end


% =========================================================================
function c = sample_pixel_corr(Mv)
%SAMPLE_PIXEL_CORR  Mean absolute inter-pixel correlation on a pixel sample.
%   Mv is [validPix x frames].  Samples up to 64 non-flat pixels and returns
%   the mean magnitude of the off-diagonal correlations (1 => collapsed).
    nP = size(Mv, 1);
    ns = min(64, nP);
    idx = round(linspace(1, nP, ns));
    S = Mv(idx, :);
    S = S(var(S, 0, 2) > 0, :);                % drop flat pixels
    if size(S, 1) < 3
        c = 0; return;                          % too few to judge
    end
    Rc = corrcoef(S');                          % pixels x pixels
    mask = triu(true(size(Rc)), 1);
    c = mean(abs(Rc(mask)), 'omitnan');
end
