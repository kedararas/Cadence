function [mask, info] = create_snr_mask(snr, floorSNR, doClean)
%CREATE_SNR_MASK  Tissue mask from an SNR map by an absolute or adaptive SNR floor.
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
%               Pass [] or 'adaptive' for a PER-RECORDING adaptive floor:
%                   floorSNR = max(1.5, 0.5 * median(snr(snr > 1.5)))
%               i.e. keep pixels at least half as bright as typical tissue,
%               never below 1.5. A fixed absolute floor does not transfer
%               across preps: on a dim heart (median tissue SNR ~3) a floor
%               of 5 keeps only the brightest corner and distorts wavefront
%               geometry (origin, planarity), while on a bright heart a floor
%               of 1.5 admits border noise that corrupts the LAT plane fit
%               and CV pool. The adaptive rule tracks each recording's own
%               brightness and reproduces careful hand-drawn masks.
%     doClean   true/false (default false).  Morphological cleanup: fill holes
%               and drop isolated specks (needs Image Processing Toolbox).
%
%   OUTPUTS
%     mask   [rows x cols] logical, true = tissue (SNR >= floorSNR, finite)
%     info   struct: .threshold (resolved value), .mode ('absolute'|'adaptive'),
%            .tissue_fraction
%
%   An absolute floor is the robust default for SNR masking: unlike Otsu /
%   triangle it makes no bimodality assumption, so it does not carve up the
%   tissue distribution when the field of view is mostly tissue. The adaptive
%   mode keeps that robustness (median-based, no bimodality assumption) while
%   scaling with per-recording illumination.

    if nargin < 2, floorSNR = 3;     end
    if nargin < 3 || isempty(doClean),  doClean  = false; end

    MIN_FLOOR = 1.5;
    adaptive = isempty(floorSNR) || ...
        ((ischar(floorSNR) || isstring(floorSNR)) && strcmpi(floorSNR, 'adaptive'));
    if adaptive
        base = isfinite(snr) & (snr > MIN_FLOOR);
        if any(base(:))
            floorSNR = max(MIN_FLOOR, 0.5 * median(snr(base)));
        else
            floorSNR = MIN_FLOOR;
        end
    end

    mask = isfinite(snr) & (snr >= floorSNR);   % NaN/Inf excluded

    if doClean
        mask = imfill(mask, 'holes');
        mask = bwareaopen(mask, max(8, round(5e-4 * numel(mask))));  % drop specks
    end

    if adaptive, mode = 'adaptive'; else, mode = 'absolute'; end
    info = struct('threshold', floorSNR, 'mode', mode, ...
                  'tissue_fraction', sum(mask(:)) / max(sum(isfinite(snr(:))), 1));

    fprintf('create_snr_mask (%s): SNR >= %.2f  ->  tissue = %.0f%% of finite pixels\n', ...
            mode, floorSNR, 100*info.tissue_fraction);
end
