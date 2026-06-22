function [start_frames, end_frames] = extract_beat_frames(pacing, varargin)
% EXTRACT_BEAT_FRAMES  Extract start and end frame indices for each beat
%   from a pacing stimulus vector.
%
%   [start_frames, end_frames] = extract_beat_frames(pacing)
%   [start_frames, end_frames] = extract_beat_frames(pacing, Name, Value)
%
%   Inputs
%     pacing        Binary vector with 1 at each beat onset (e.g. from
%                   auto_detect_peaks), OR an analog pacing trace with
%                   voltage pulses.  Length must equal the number of
%                   frames in the recording.
%
%   Name-Value options
%     'PreWindow'   Frames to include BEFORE each beat onset (default 0).
%                   start_frames(i) = onset(i) - PreWindow.
%                   Useful for capturing a pre-stimulus baseline.
%
%     'FixedLength' true | false (default false).
%                   false : each beat spans onset(i) to onset(i+1)-1;
%                           the last beat spans onset(end) to T.
%                   true  : every beat (including the last) uses the
%                           median inter-beat interval as the window
%                           length, so all windows are the same size.
%                           Useful for ensemble averaging.
%
%   Outputs
%     start_frames  [n_beats x 1] first frame index of each beat window.
%     end_frames    [n_beats x 1] last  frame index of each beat window.
%
%   Notes
%     • Beats whose window would extend before frame 1 or after frame T
%       are flagged with a warning and excluded from the output.
%     • start_frames and end_frames are 1-based frame indices that can
%       be used directly for indexing: data(:, :, start:end).
%
%   Example
%     pacing = auto_detect_peaks(data);
%     [sf, ef] = extract_beat_frames(pacing, 'PreWindow', 10, 'FixedLength', true);
%     for k = 1:numel(sf)
%         beat_k = data(:, :, sf(k):ef(k));
%     end

    % ── Parse inputs ──────────────────────────────────────────────────────
    p = inputParser;
    addRequired(p,  'pacing',      @isnumeric);
    addParameter(p, 'PreWindow',   0,     @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'FixedLength', false, @(x) islogical(x) || isnumeric(x));
    parse(p, pacing, varargin{:});

    pre_win      = round(p.Results.PreWindow);
    fixed_length = logical(p.Results.FixedLength);
    T            = numel(pacing);

    % ── Find beat onsets ──────────────────────────────────────────────────
    if all(pacing == 0 | pacing == 1)
        % Binary vector — onsets are the 1-valued frames
        onsets = find(pacing(:));
    else
        % Analog trace — find stimulus pulses
        min_sep = round(T / 10);
        [~, onsets] = findpeaks(double(pacing(:)), ...
            'MinPeakHeight',   0.5 * max(pacing), ...
            'MinPeakDistance', min_sep);
    end

    if isempty(onsets)
        warning('extract_beat_frames:noOnsets', ...
            'No beat onsets found in pacing vector. Returning empty arrays.');
        start_frames = zeros(0, 1);
        end_frames   = zeros(0, 1);
        return;
    end

    n_beats = numel(onsets);

    % ── Compute inter-beat intervals ──────────────────────────────────────
    ibi        = diff(onsets);          % length n_beats-1
    median_ibi = median(ibi);

    % ── Build start / end frames ──────────────────────────────────────────
    sf = onsets - pre_win;              % shift back by pre-stimulus window

    if fixed_length
        % Every beat gets the same window length (median IBI)
        beat_len = round(median_ibi);
        ef = sf + beat_len - 1;
    else
        % Each beat spans to just before the next onset
        % Last beat: extend by median IBI (best estimate for unknown next onset)
        ef        = [onsets(2:end) - 1; onsets(end) + round(median_ibi) - 1];
        % Shift end frames back if we shifted start frames back
        % (the window length stays onset-to-onset regardless of pre_win)
        % ef already points to onset(i+1)-1, which is correct —
        % pre_win shifts start back without changing where the beat ends.
    end

    % ── Clip to recording bounds and warn about out-of-range windows ──────
    in_range = sf >= 1 & ef <= T;

    if any(~in_range)
        n_dropped = sum(~in_range);
        warning('extract_beat_frames:windowOutOfRange', ...
            '%d beat(s) dropped: window extends outside recording [1, %d].', ...
            n_dropped, T);
    end

    start_frames = sf(in_range);
    end_frames   = ef(in_range);

    % ── Summary ───────────────────────────────────────────────────────────
    fprintf('extract_beat_frames: %d beats detected, %d valid.\n', ...
        n_beats, numel(start_frames));
    fprintf('  Median IBI : %d frames\n', round(median_ibi));
    fprintf('  Window     : %d frames (%s)\n', ...
        round(median_ibi) + pre_win, ...
        ternary(fixed_length, 'fixed length', 'onset-to-onset'));
end


% ── Local helper ──────────────────────────────────────────────────────────
function s = ternary(cond, a, b)
    if cond; s = a; else; s = b; end
end
