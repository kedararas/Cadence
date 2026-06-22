function cmos_all_data = extract_ps_dynamics(cmos_all_data, min_cycles, threshold_mode, target_fdr)
%EXTRACT_PS_DYNAMICS  Track phase singularities and flag valid (rotor) PS.
%
%   cmos_all_data = extract_ps_dynamics(cmos_all_data)
%   cmos_all_data = extract_ps_dynamics(cmos_all_data, min_cycles)
%   cmos_all_data = extract_ps_dynamics(cmos_all_data, min_cycles, threshold_mode, target_fdr)
%
%   min_cycles      Minimum number of full rotation periods a PS must persist
%                   to be flagged as a valid rotor (default 1.5).  Used in
%                   'cycles' mode, and as the fallback in 'fdr' mode.  At ~1.5
%                   cycles the estimated false-discovery rate (chance-linked
%                   transient PS) is ~10-15%; at 1.0 it is ~50%; at 2.0 it is
%                   high confidence but very few PS qualify.
%
%   threshold_mode  'cycles' (default) | 'fdr'.
%                   'cycles' — validity threshold = min_cycles rotation periods
%                              (a fixed, frequency-scaled cycle count).
%                   'fdr'    — validity threshold is chosen from the data so the
%                              estimated false-discovery rate of the surviving
%                              PS falls at/below target_fdr.  Self-calibrates to
%                              each recording's noise level; falls back to the
%                              min_cycles threshold if it cannot be fit.  See
%                              estimate_fdr_threshold.
%
%   target_fdr      Target false-discovery rate for 'fdr' mode (default 0.10).
%
%   Validate either choice against expert visual scoring before relying on it
%   across species/recordings.
if nargin < 2 || isempty(min_cycles)
    min_cycles = 1.5;
end
if nargin < 3 || isempty(threshold_mode)
    threshold_mode = 'cycles';
end
if nargin < 4 || isempty(target_fdr)
    target_fdr = 0.10;
end

% Detect and annotate PS — two thresholds are computed inside count_ps:
%   dup_threshold:      merges nearby detections of the same PS (cluster radius)
%   wf_assoc_threshold: how close a PS must be to a wavefront tip to be tagged
% Pass pixel_size_mm if known for physically-grounded thresholds; omit for
% geometry-derived fallbacks.
pixel_size_mm = cmos_all_data.pixel_size;   % set to e.g. 0.1 if known (mm per pixel)
[ps, ps_count] = count_ps(cmos_all_data.ps_data, cmos_all_data.wavefronts, pixel_size_mm);
cmos_all_data.ps       = ps;
cmos_all_data.ps_count = ps_count;

% Minimum PS lifespan threshold in frames.
% A PS must persist for at least min_cycles full rotation periods to be valid
% (min_cycles is the optional function argument above; default 1.5).
% Requires frame_rate (frames/sec) and a rotor-frequency estimate (Hz).
frame_rate = cmos_all_data.frame_rate;   % frames/sec — must be present in cmos_all_data

% Rotor frequency for the rotation period.  See estimate_rotor_frequency:
% the DF map is bimodal, so a high percentile (not the median) is used so the
% period reflects the fast, rotor-bearing tissue rather than slow far-field
% pixels.  Shared with the wavefront-dynamics path for a single definition.
df_percentile = 75;
df_rotor      = estimate_rotor_frequency(cmos_all_data.df_map, df_percentile);  % Hz
cycle_frames  = frame_rate / df_rotor;                    % frames per rotation
cycles_threshold = round(min_cycles * cycle_frames);      % fixed-cycle threshold
% The validity threshold is applied AFTER tracking (see end of function).
% In 'fdr' mode cycles_threshold is the fallback if the noise fit fails.

% Look-ahead window: allow PS to be missing for up to 25% of one rotation period.
% Expressed in frames so it scales correctly with frame rate.
% Tuning guide: increase missing_fraction if rotors are frequently split into
% short tracks; decrease if unrelated PS are being incorrectly linked.
missing_fraction = 0.25;
lookahead_frames = round(missing_fraction * frame_rate / df_rotor);

% PS tracking threshold — maximum deviation from predicted position (pixels).
% With velocity extrapolation, this covers detection noise only (not drift),
% so it can be tight. During look-ahead it widens with sqrt(frames_ahead).
% Tuning guide: if valid PS are being terminated too early, increase by 1-2 px;
% if unrelated PS are being linked, decrease.
num_ps_rows = size(cmos_all_data.ps_data, 1);
num_ps_cols = size(cmos_all_data.ps_data, 2);
if ~isempty(pixel_size_mm) && pixel_size_mm > 0
    % Physical: allow up to 1 mm deviation from predicted position
    ps_track_threshold = round(1.0 / pixel_size_mm);
else
    % Geometry fallback: ~1.5% of frame diagonal
    ps_track_threshold = max(3, round(0.015 * sqrt(num_ps_rows^2 + num_ps_cols^2)));
end
fprintf(['extract_ps_dynamics: df_rotor(p%d)=%.2f Hz (%.0f frames/cycle), ' ...
    'ps_track_threshold=%d px, lookahead=%d frames, mode=%s\n'], ...
    df_percentile, df_rotor, cycle_frames, ps_track_threshold, lookahead_frames, ...
    lower(threshold_mode));

% Build PS object matrix from row 1 only — all cleaned PS with wavefront metadata.
% ps(1,:) = [x, y, charge, wf_dist, wf_index] per frame (5 columns).
% Passing the full 3×num_frames cell would cause numel to return 3× the frame
% count, iterating over raw and cleaned rows incorrectly.
[ps_matrix, valid_ps] = create_ps_objects(cmos_all_data.ps(1,:));

% PSDynamics is a handle class — no struct copy on each update call
ps_dynamics = PSDynamics(valid_ps);

num_frames = size(ps_matrix, 2);

for i = 1:num_frames
    ps_list = ps_matrix{i};
    if isempty(ps_list)
        continue;
    end

    num_ps = numel(ps_list);

    for j = 1:num_ps  % iterate through every ps in a frame
        ps = ps_list{j};

        % Set birthday if needed; skip already-tracked singularities
        if isempty(ps.birthday)
            ps.birthday     = [ps.frame, ps.index];
            ps.lifespan     = [i, 0, i];
            ps.distance     = 0;
            ps.displacement = 0;
        else
            continue;  % no double dipping
        end

        current_ps_list = ps_list;  % start comparison from current frame
        k = i + 1;

        % Pre-allocate path buffers — O(1) per-frame append, one final slice copy.
        % Growing by concatenation ([path; new_row]) is O(N²) total for N-frame PS.
        max_life = num_frames - i + 1;
        path_buf = zeros(max_life, 4);  % [frame, index, x, y]
        wf_buf   = zeros(max_life, 2);  % [frame, wf_index]
        path_len = 1;
        path_buf(1,:) = [i, j, ps.location];
        wf_buf(1,:)   = ps.wavefront;

        ps_terminated = false;
        while k <= num_frames
            if ~isempty(ps_matrix{k})
                next_ps_list = ps_matrix{k};

                % Predict where this PS should be in frame k using recent velocity
                expected_pos = predict_ps_position(path_buf, path_len, k);
                [ps_status, ps_cand] = choose_ps_cand(ps, current_ps_list, next_ps_list, ...
                    ps_track_threshold, expected_pos);

                % Look-ahead: if not found, search ahead up to lookahead_frames.
                % Threshold widens with sqrt(frames_ahead) to account for growing
                % extrapolation uncertainty, capped at 3x the base threshold.
                if strcmp(ps_status, 'EXPIRED') && k+1 <= num_frames
                    end_window = min(k + lookahead_frames, num_frames);
                    for m = k+1:end_window
                        frames_ahead     = m - path_buf(path_len, 1);
                        scaled_threshold = min(ps_track_threshold * sqrt(frames_ahead), ...
                                               ps_track_threshold * 3);
                        expected_pos = predict_ps_position(path_buf, path_len, m);
                        [ps_status, ps_cand] = choose_ps_cand(ps, current_ps_list, ps_matrix{m}, ...
                            scaled_threshold, expected_pos);
                        if strcmp(ps_status, 'ALIVE')
                            k = m;
                            break;
                        end
                    end
                end

                switch ps_status
                    case 'ALIVE'
                        ps.lifespan(1,3) = k;
                        ps.lifespan(1,2) = k - ps.lifespan(1,1);
                        path_len = path_len + 1;
                        path_buf(path_len,:) = [ps_cand.frame, ps_cand.index, ps_cand.location];
                        wf_buf(path_len,:)   = ps_cand.wavefront;
                        ps.distance = ps.distance + hypot( ...
                            ps_cand.location(1) - ps.location(1), ...
                            ps_cand.location(2) - ps.location(2));
                        ps.displacement = hypot( ...
                            path_buf(path_len,3) - path_buf(1,3), ...
                            path_buf(path_len,4) - path_buf(1,4));

                        current_ps_list      = next_ps_list;
                        ps_cand.lifespan     = ps.lifespan;
                        ps_cand.distance     = ps.distance;
                        ps_cand.displacement = ps.displacement;
                        ps_cand.birthday     = ps.birthday;
                        ps = ps_cand;
                        k = k + 1;

                    case 'EXPIRED'
                        ps.path    = path_buf(1:path_len, :);
                        ps.wf_path = wf_buf(1:path_len, :);
                        ps_matrix  = backfill_ps_path(ps_matrix, ps);
                        ps_dynamics = update_ps_database(ps_dynamics, ps);
                        ps_terminated = true;
                        break;
                end

            else
                % Empty frame — PS expired due to gap
                ps.path    = path_buf(1:path_len, :);
                ps.wf_path = wf_buf(1:path_len, :);
                ps_matrix  = backfill_ps_path(ps_matrix, ps);
                ps_dynamics = update_ps_database(ps_dynamics, ps);
                ps_terminated = true;
                break;
            end
        end

        % PS alive through the last frame — record it
        if ~ps_terminated
            ps.path    = path_buf(1:path_len, :);
            ps.wf_path = wf_buf(1:path_len, :);
            ps_matrix  = backfill_ps_path(ps_matrix, ps);
            ps_dynamics = update_ps_database(ps_dynamics, ps);
        end
    end
end

% ── Apply the validity threshold (after tracking, so the full lifespan ──────
% distribution is available for the data-driven FDR mode). ──────────────────
switch lower(threshold_mode)
    case 'fdr'
        [valid_ps_threshold, fdr_info] = estimate_fdr_threshold( ...
            ps_dynamics.ps_info(:,3), cycle_frames, target_fdr, cycles_threshold);
        fprintf(['extract_ps_dynamics: [fdr mode] valid_ps_threshold=%d frames ' ...
            '(%.2f cycles) — %s\n'], ...
            valid_ps_threshold, valid_ps_threshold / cycle_frames, fdr_info.reason);
    otherwise   % 'cycles'
        valid_ps_threshold = cycles_threshold;
        fprintf(['extract_ps_dynamics: [cycles mode] valid_ps_threshold=%d frames ' ...
            '(%.2f cycles)\n'], valid_ps_threshold, min_cycles);
end

ps_dynamics = finalize_ps_validity(ps_dynamics, valid_ps_threshold);

cmos_all_data.ps_dynamics = ps_dynamics;

end  % end of main function


%% Local helper — predict PS position in frame target_frame via linear extrapolation
function expected_pos = predict_ps_position(path_buf, path_len, target_frame)
% Uses the last 3 known positions to estimate velocity (mean displacement/frame),
% then projects forward to target_frame. Falls back to last known position
% if fewer than 2 path entries exist (no velocity estimate possible).
if path_len >= 3
    recent   = path_buf(path_len-2:path_len, 3:4);  % last 3 [x, y] positions
    velocity = mean(diff(recent), 1);                % mean displacement per frame
elseif path_len >= 2
    velocity = path_buf(path_len,3:4) - path_buf(path_len-1,3:4);
else
    velocity = [0, 0];                               % stationary assumption
end
frames_ahead = target_frame - path_buf(path_len, 1);
expected_pos = path_buf(path_len, 3:4) + velocity * frames_ahead;
end


%% Local helper — writes terminal PS state back to every frame in ps.path
function ps_matrix = backfill_ps_path(ps_matrix, ps)
if isempty(ps.path)
    return;
end
for p = 1:size(ps.path, 1)
    frame = ps.path(p, 1);
    index = ps.path(p, 2);
    entry = ps_matrix{frame}{index};
    entry.lifespan     = ps.lifespan;
    entry.path         = ps.path;
    entry.distance     = ps.distance;
    entry.displacement = ps.displacement;
    entry.wf_path      = ps.wf_path;
    entry.birthday     = ps.birthday;
    ps_matrix{frame}{index} = entry;
end
end
