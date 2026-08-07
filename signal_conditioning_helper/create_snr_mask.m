function [mask, info] = create_snr_mask(snr, floorSNR, doClean)
%CREATE_SNR_MASK  Tissue mask from an SNR map by an absolute SNR floor.
%
%   mask = create_snr_mask(snr)
%   mask = create_snr_mask(snr, floorSNR)
%   [mask, info] = create_snr_mask(snr, floorSNR, doClean)
%
%   INPUTS
%     snr       [rows x cols] per-pixel SNR map (e.g. from extract_snr_mask)
%     floorSNR  keep pixels with SNR >= floorSNR (default 3).  ~3-5 is the
%               conventional "signal reliably above noise" range; use ~3 for a
%               permissive processing mask (SVD / pooling) and ~5 for a
%               stricter analysis mask (APD / feature extraction).
%     doClean   true/false (default false).  Morphological cleanup: fill holes
%               and drop isolated specks (needs Image Processing Toolbox).
%
%   OUTPUTS
%     mask   [rows x cols] logical, true = tissue (SNR >= floorSNR, finite)
%     info   struct: .threshold, .tissue_fraction
%
%   An absolute floor is the robust default for SNR masking: unlike Otsu /
%   triangle it makes no bimodality assumption, so it does not carve up the
%   tissue distribution when the field of view is mostly tissue.

    if nargin < 2 || isempty(floorSNR), floorSNR = 3;     end
    if nargin < 3 || isempty(doClean),  doClean  = false; end

    mask = isfinite(snr) & (snr >= floorSNR);   % NaN/Inf excluded

    if doClean
        mask = imfill(mask, 'holes');
        mask = bwareaopen(mask, max(8, round(5e-4 * numel(mask))));  % drop specks
    end

    info = struct('threshold', floorSNR, ...
                  'tissue_fraction', sum(mask(:)) / max(sum(isfinite(snr(:))), 1));

    fprintf('create_snr_mask: SNR >= %g  ->  tissue = %.0f%% of finite pixels\n', ...
            floorSNR, 100*info.tissue_fraction);
end
