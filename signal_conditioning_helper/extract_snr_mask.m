function snr_mask = extract_snr_mask(cmos_data)
%EXTRACT_SNR_MASK  Compute per-pixel SNR as signal amplitude / noise std dev.
%
%   snr_mask = extract_snr_mask(cmos_data)
%
%   Formula
%     SNR(i,j) = amplitude(i,j) / noise_std(i,j)
%
%     amplitude = max(smooth) - min(smooth)   peak-to-peak of the true OAP
%     noise_std = std(raw - smooth)           std dev of the residual noise
%
%   Both terms have units of fluorescence counts, so the result is
%   dimensionless and scale-invariant — the standard definition used in
%   the cardiac optical mapping literature (Laughner et al. and others).
%
%   The smoothed trace is estimated with a centred moving average (movmean)
%   whose window is 1 % of the total frame count (minimum 3 frames).
%   movmean operates natively on the 3-D array along the time dimension,
%   so no pixel loop is needed.  The residual (raw − smooth) captures the
%   high-frequency noise; its std dev is the noise floor estimate.
%   Pixels with zero noise std (dead / uniform background) are left as 0.
%
%   Input
%     cmos_data   [rows × cols × frames]  raw optical-mapping data
%   Output
%     snr_mask    [rows × cols]            SNR map (amplitude / noise std dev)

[num_rows, num_cols, T] = size(cmos_data);

% ── Step 1: Smooth along the time dimension (no pixel loop) ──────────────────
% Window = 1 % of frame count, minimum 3 frames.
% movmean(A, k, 3) applies a centred k-frame moving average along dim 3
% to the entire rows×cols×T array in one vectorized call.
win         = max(3, round(0.01 * T));
smooth_data = movmean(double(cmos_data), win, 3);

% ── Step 2: Compute SNR map (fully vectorized) ────────────────────────────────
residual  = double(cmos_data) - smooth_data;                   % rows × cols × T
amplitude = max(smooth_data, [], 3) - min(smooth_data, [], 3); % rows × cols
noise_std = std(residual, 0, 3);                                % rows × cols

% Initialise to zero; only fill pixels that have non-zero noise std.
% Pixels with noise_std == 0 are dead / uniform background — SNR left as 0.
snr_mask        = zeros(num_rows, num_cols);
valid           = noise_std > 0;
snr_mask(valid) = amplitude(valid) ./ noise_std(valid);

end
