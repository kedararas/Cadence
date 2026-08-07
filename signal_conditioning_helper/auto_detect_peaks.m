function pacing = auto_detect_peaks(data)
% AUTO_DETECT_PEAKS  Estimate beat onset times from a 3-D optical-mapping
%   stack when no pacing stimulus vector (analog1) is available.
%
%   pacing = auto_detect_peaks(data)
%
%   Returns a binary vector with 1 at each detected beat onset frame.
%
%   Input
%     data   [rows x cols x T] optical-mapping array (any numeric class)
%
%   Strategy
%     1. Collapse to a spatial-mean 1-D signal (one scalar per frame).
%     2. Smooth with a Gaussian (2 % of recording, min 5 frames).
%     3. Normalise to [0, 1].
%     4. Compute the first derivative via diff.
%     5. Find peaks in the derivative — the frame of maximum upstroke
%        rate is earlier than the signal peak and marks beat onset.
%     6. Enforce a minimum inter-beat separation (T/200 frames) to
%        suppress spurious peaks on notches and the falling edge.
%
%   NOTE: The definitive fix for stimulus alignment is to pass the
%   analog1 pacing vector from GSDconverter_BVW directly to
%   ensembleAverageFull rather than calling this function.

    [~, ~, T] = size(data);

    % ── Step 1: Spatial mean ───────────────────────────────────────────────
    sig = reshape(mean(reshape(double(data), [], T), 1), T, 1);

    % ── Step 2: Smooth (Gaussian, 2 % of recording, min 5 frames) ─────────
    win = max(5, round(T * 0.02));
    sig = smoothdata(sig, 'gaussian', win);

    % ── Step 3: Normalise to [0, 1] ───────────────────────────────────────
    sig_min = min(sig);
    sig_max = max(sig);
    if sig_max > sig_min
        sig = (sig - sig_min) / (sig_max - sig_min);
    else
        warning('auto_detect_peaks:flatSignal', 'Signal has no range — returning empty pacing.');
        pacing = zeros(T, 1);
        return;
    end

    % ── Step 4: Derivative via diff ───────────────────────────────────────
    dsig = diff(sig);   % length T-1

    % ── Step 4b: Handle inverted signals ──────────────────────────────────
    % Some dye/camera combinations produce a downward deflection at beat
    % onset (fluorescence decreases on depolarisation). In that case the
    % derivative at onset is a large NEGATIVE value, not a large positive
    % one, and findpeaks would find the wrong phase (recovery upswing).
    % Detection: if the largest-magnitude derivative value is negative,
    % the signal is inverted — negate so findpeaks always sees upstrokes.
    if abs(min(dsig)) > abs(max(dsig))
        dsig = -dsig;
    end

    % ── Step 5: Find peaks in dF/dt ───────────────────────────────────────
    min_sep_frames  = round(T / 200);
    min_peak_height = 0.10 * max(dsig);

    [~, locs] = findpeaks(dsig, ...
        'MinPeakHeight',   min_peak_height, ...
        'MinPeakDistance', min_sep_frames);

    if isempty(locs)
        warning('auto_detect_peaks:noPeaksFound', ...
            ['Auto peak detection found no beats.\n' ...
             'Consider passing the analog1 pacing vector explicitly.']);
        pacing = zeros(T, 1);
        return;
    end

    % ── Step 6: Remove erroneous peaks using median IBI ───────────────────
    % A spurious extra detection creates a short IBI followed immediately
    % by another short IBI where one normal interval should be.
    % Strategy: iteratively find the shortest IBI; if it is less than 50 %
    % of the current median, the pair shares a double-detection — drop
    % whichever peak has the smaller dsig value (weaker upstroke).
    % Repeat until no short intervals remain.
    if numel(locs) >= 3
        while true
            ibi     = diff(locs);
            med_ibi = median(ibi);
            [min_ibi, idx] = min(ibi);          % shortest interval

            if min_ibi >= 0.5 * med_ibi         % all intervals look normal
                break;
            end

            % The pair is locs(idx) and locs(idx+1)
            % Remove whichever has the smaller derivative peak
            peak_a = dsig(locs(idx));
            peak_b = dsig(locs(idx + 1));
            if peak_a >= peak_b
                locs(idx + 1) = [];             % drop the weaker second peak
            else
                locs(idx)     = [];             % drop the weaker first peak
            end
        end
    end

    pacing       = zeros(T, 1);
    pacing(locs) = 1;
end
