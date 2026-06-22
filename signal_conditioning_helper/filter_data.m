function [filt_data, Wp] = filter_data(data,Fs,Wp,method,sgOrder)
%FILTER_DATA  Zero-phase temporal filtering of an optical-mapping stack.
%
%   filt_data = filter_data(data, Fs, Wp)                 % Butterworth (default)
%   filt_data = filter_data(data, Fs, 'auto')             % auto-pick cutoff
%   filt_data = filter_data(data, Fs, Wp, 'sgolay')       % Savitzky-Golay
%   filt_data = filter_data(data, Fs, Wp, 'sgolay', poly) % SG, polynomial order
%   [filt_data, Wp] = filter_data(...)                    % also return Wp used
%
% INPUTS
% data    = cmos data [rows x cols x frames]
% Fs      = sampling frequency (Hz)
% Wp      = approximate low-pass cutoff frequency (Hz), OR 'auto' / [] to
%           estimate it from the data (see AUTO CUTOFF below)
% method  = 'butter' (default) | 'sgolay'
% sgOrder = Savitzky-Golay polynomial order (default 3; 'sgolay' only)
%
% AUTO CUTOFF
% With Wp = 'auto' (or []), the cutoff is estimated per recording so one call
% adapts across species/heart-rates.  A high-SNR global trace (spatial mean of
% the highest-variance pixels) is built, the beat period is found by
% autocorrelation, the beats are folded into a clean averaged beat, and the
% upstroke 10-90% rise time tr is measured.  The signal bandwidth is then
% ~0.35/tr; the cutoff is set to 1.5x that, CLAMPED to [40, 80] Hz (a range
% validated to preserve the upstroke across rat/rabbit/pig/human while
% removing out-of-band noise).  If the period/rise cannot be found it falls
% back to 50 Hz.  NOTE: on noisy raw stacks the estimate is approximate and
% typically lands near 40-50 Hz; pass an explicit Wp when you need exact
% control (e.g. a wider band to preserve a very fast upstroke).
%
% OUTPUT
% filt_data = filtered data, same size as input
%
% METHOD
% 'butter'  4th-order Butterworth IIR low-pass in second-order-section (SOS)
%           form, applied with filtfilt (forward + backward) for zero phase
%           (effective 8th-order roll-off, no group-delay distortion).
%
% 'sgolay'  Savitzky-Golay smoothing: a least-squares polynomial fit over a
%           sliding window.  It removes broadband noise while preserving peak
%           height and the up/down-stroke shape far better than a Butterworth
%           low-pass of equivalent smoothing, which matters for APD, dV/dt-max
%           and activation-time metrics.  The window length is set from Wp so
%           it still acts as an approximate cutoff: frame ~ Fs/Wp samples
%           (forced odd, and > sgOrder).  Zero-phase by construction.

if nargin < 4 || isempty(method),  method  = 'butter'; end
if nargin < 5 || isempty(sgOrder), sgOrder = 3;        end

%% Resolve the cutoff (explicit, or estimated when 'auto' / [])
if isempty(Wp) || (ischar(Wp) || isstring(Wp)) && strcmpi(Wp, 'auto')
    Wp = estimate_cutoff(data, Fs);
    fprintf('filter_data: auto cutoff Wp = %.0f Hz\n', Wp);
elseif ischar(Wp) || isstring(Wp)
    error('filter_data:badWp', 'Wp must be a number or ''auto''.');
end

[R, C, T] = size(data);
% Reshape to [T x pixels] so the 1-D filters operate column-wise on every
% pixel at once (pre-R2020b filtfilt does not support N-D arrays).
tmp = reshape(permute(data, [3,1,2]), T, R*C);

switch lower(method)
    case 'butter'
        %% Butterworth low-pass (SOS form), cached across calls
        persistent sos_cache g_cache Fs_cache Wp_cache
        if isempty(sos_cache) || Fs ~= Fs_cache || Wp ~= Wp_cache
            Wn  = Wp / (Fs / 2);                 % normalised cutoff (0..1)
            Wn  = min(Wn, 0.9999);               % guard against cutoff >= Nyquist
            [z,p,k] = butter(4, Wn, 'low');      % ZPK avoids TF coefficient errors
            [sos_cache, g_cache] = zp2sos(z, p, k);
            Fs_cache = Fs;
            Wp_cache = Wp;
        end
        tmp = filtfilt(sos_cache, g_cache, tmp);

    case 'sgolay'
        %% Savitzky-Golay smoothing
        % Window spanning ~1/Wp seconds, forced odd and large enough for the
        % polynomial order.  sgolayfilt is zero-phase (symmetric window).
        frame = round(Fs / Wp);
        if mod(frame, 2) == 0, frame = frame + 1; end   % must be odd
        frame = max(frame, sgOrder + 2 + mod(sgOrder,2)); % > order, and odd
        if mod(frame, 2) == 0, frame = frame + 1; end
        frame = min(frame, T - mod(T+1,2));             % cannot exceed length
        tmp = sgolayfilt(tmp, sgOrder, frame);

    otherwise
        error('filter_data:badMethod', ...
              'Unknown method "%s". Use ''butter'' or ''sgolay''.', method);
end

filt_data = permute(reshape(tmp, T, R, C), [2,3,1]);
end


% =========================================================================
function Wp = estimate_cutoff(data, Fs)
%ESTIMATE_CUTOFF  Per-recording low-pass cutoff from the upstroke kinetics.
%   Bandwidth ~ 0.35 / (10-90% upstroke rise time), measured on a clean
%   ensemble-averaged beat of a high-SNR pooled trace.  Clamped to a safe
%   band; falls back to 50 Hz when the beat/upstroke cannot be resolved.

    WP_MIN = 40; WP_MAX = 80; WP_FALLBACK = 50; MARGIN = 1.5;

    [R, C, T] = size(data);
    M = reshape(permute(double(data), [3 1 2]), T, R*C);   % T x pixels

    % High-SNR pooled trace: spatial mean of the top-quartile-variance pixels
    sd  = std(M, 0, 1);
    roi = sd >= prctile(sd(isfinite(sd)), 75);
    g   = mean(M(:, roi & isfinite(sd)), 2, 'omitnan');    % T x 1

    % Remove bleach/drift (5th-order polynomial)
    tt = (1:T)';
    g  = g - polyval(polyfit(tt, g, 5), tt);

    % Beat period: STRONGEST autocorrelation peak in a physiological lag range
    g0 = g - mean(g);
    ac = xcorr(g0, 'normalized');
    ac = ac(T:end);                                        % lags 0..T-1
    loLag = 50; hiLag = min(1500, floor(T/3));
    seg = ac(loLag+1 : hiLag+1);
    [pks, locs] = findpeaks(seg);
    strong = pks > 0.2;
    if ~any(strong), Wp = WP_FALLBACK; return; end
    locs = locs(strong); pks = pks(strong);
    [~, im] = max(pks);
    CL = locs(im) + loLag;                                 % cycle length (samples)

    % Fold beats into one clean averaged beat
    nb   = floor(T / CL);
    if nb < 2, Wp = WP_FALLBACK; return; end
    fold = mean(reshape(g(1:nb*CL), CL, nb), 2);
    fold = (fold - min(fold)) / (max(fold) - min(fold) + eps);
    if mean(fold > 0.5) > 0.5, fold = 1 - fold; end        % ensure AP upward

    % 10-90% upstroke rise time around the steepest point
    [~, u] = max(diff(fold));
    a = max(1, u-10); b = min(CL, u+30);
    s = fold(a:b); s = (s - min(s)) / (max(s) - min(s) + eps);
    i10 = find(s >= 0.1, 1, 'first');
    i90 = find(s >= 0.9, 1, 'first');
    if isempty(i10) || isempty(i90) || i90 <= i10
        Wp = WP_FALLBACK; return;
    end
    tr = (i90 - i10) / Fs;                                 % rise time (s)
    Wp = min(max(0.35/tr * MARGIN, WP_MIN), WP_MAX);
end
