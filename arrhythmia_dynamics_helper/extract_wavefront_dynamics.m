function wavefront_dynamics = extract_wavefront_dynamics(cmos_all_data,threshold,debug_mode)
%EXTRACT_WAVEFRONT_DYNAMICS  Track activation wavefronts across frames.
%
%   wavefront_dynamics = extract_wavefront_dynamics(cmos_all_data, threshold, debug_mode)
%
%   Top-level wavefront-tracking entry point (called from the Feature
%   Extraction module). Builds Wavefront objects from the per-frame segmented
%   wavefront list, then frame-by-frame links each wavefront to its best match
%   in the next frame (via get_wf_neighbors / choose_wf_cand) and records the
%   resulting tracks in a WavefrontDynamics handle object. The dominant
%   frequency and frame rate carried on cmos_all_data seed the dynamics object.
%
%   Inputs
%     cmos_all_data  struct with fields wavefronts (2 x num_frames cell;
%                    row 1 = tracking-quality segments), wf_count, df_map,
%                    frame_rate.
%     threshold      vector of tracking parameters; threshold(1,7) is the
%                    neighbour / matching distance threshold (pixels).
%     debug_mode     flag enabling diagnostic plots / console output.
%
%   Output
%     wavefront_dynamics  WavefrontDynamics object holding all tracked
%                         wavefront paths and summary statistics.


try

% WavefrontDynamics is a handle class — no struct copy overhead as it grows
wavefront_dynamics = WavefrontDynamics( ...
    cmos_all_data.wf_count, ...
    median(cmos_all_data.df_map, 'all', 'omitnan'), ...
    cmos_all_data.frame_rate);

% Wavefront Tracking
% cmos_all_data.wavefronts is 2×num_frames: row 1 = filtered newsegWave (for
% tracking), row 2 = all non-empty segments (for PS tagging in count_ps).
% Pass row 1 only — numel on the full 2×N cell would return 2*num_frames and
% linear indexing would interleave tracking-quality and raw-segment frames.
wavefronts = create_wavefront_objects(cmos_all_data.wavefronts(1,:), 6*threshold(1,7));
num_frames = size(wavefronts, 2);
dist_thresh = threshold(1,7); % cache to avoid repeated indexing

for i = 1:num_frames % iterate through every frame
    if mod(i, 100) == 0 || i == 1
        fprintf('Processing frame %d / %d\n', i, num_frames);
    end

    wf_list = wavefronts{1,i};
    if isempty(wf_list)
        continue;
    end

    num_wf = numel(wf_list);

    for j = 1:num_wf % iterate through every wavefront in a frame
        wf = wf_list{j,1};

        % Set birthday if needed; skip already-tracked wavefronts
        if isempty(wf.birthday)
            wf.birthday = [wf.frame, wf.index];
            wf.lifespan = [i, 0, i];
            current_neighbors = get_wf_neighbors(wf, wf_list, dist_thresh, debug_mode);
            if ~isempty(current_neighbors)
                current_neighbors(:,3) = wf.frame;
                current_neighbors(current_neighbors(:,2) == wf.index, :) = [];
                if ~isempty(current_neighbors)
                    wf.neighbors{1,1} = current_neighbors;
                end
            end
        else
            continue; % no double dipping
        end

        neighbor_count = numel(wf.neighbors); % track neighbors list length
        current_wf_list = wf_list;            % start comparison from current frame

        % Pre-allocate path and length buffers — O(1) per-frame append,
        % one final slice copy at termination. Growing by concatenation
        % ([wf.path; new_row]) is O(N²) total for an N-frame wavefront.
        max_life = num_frames - i + 1;
        path_buf = zeros(max_life, 2);  % [frame, index]
        len_buf  = zeros(max_life, 1);  % wavefront length per frame
        path_len = 1;
        path_buf(1,:) = [i, j];
        len_buf(1)    = wf.length;      % scalar from create_wavefront_objects

        wf_terminated = false;
        for k = i+1:num_frames % compare with future frames

            if ~isempty(wavefronts{1,k})
                next_wf_list = wavefronts{1,k};
                next_wf_black_list = [];

                [wf_status, wf_cand, wf_mate, wf_child] = choose_wf_cand(wf, current_wf_list, next_wf_list, next_wf_black_list, dist_thresh, debug_mode);

                switch wf_status
                    case 'ALIVE'
                        neighbors = get_wf_neighbors(wf_cand, next_wf_list, dist_thresh, debug_mode);
                        if ~isempty(neighbors)
                            neighbors(:,3) = wf_cand.frame;
                            neighbors(neighbors(:,2) == wf_cand.index, :) = [];
                        end

                        wf.lifespan(1,2) = wf.lifespan(1,2) + 1;
                        wf.lifespan(1,3) = k;
                        path_len = path_len + 1;
                        path_buf(path_len,:) = [wf_cand.frame, wf_cand.index];
                        len_buf(path_len)    = median(wf_cand.length);
                        if ~isempty(neighbors)
                            neighbor_count = neighbor_count + 1;
                            wf.neighbors{neighbor_count, 1} = neighbors;
                        end

                        current_wf_list   = next_wf_list;
                        wf_cand.lifespan  = wf.lifespan;
                        wf_cand.length    = len_buf(path_len);  % scalar: current length for choose_wf_cand
                        wf_cand.neighbors = wf.neighbors;
                        wf_cand.birthday  = wf.birthday;
                        wf = wf_cand;

                    case 'EXPIRED'
                        wf.path           = path_buf(1:path_len, :);
                        wf.length         = len_buf(1:path_len);
                        wf.cause_of_death = 'EXPIRED';
                        wavefronts        = backfill_wf_path(wavefronts, wf, struct());
                        wavefront_dynamics = update_wf_database(wavefront_dynamics, wf, wavefronts, threshold);
                        wf_terminated = true;
                        break;

                    case 'FRAGMENTED'
                        wf.path           = path_buf(1:path_len, :);
                        wf.length         = len_buf(1:path_len);
                        wf.cause_of_death = 'FRAGMENTED';
                        % Register each child's parent reference
                        for c = 1:numel(wf_child)
                            child = wf_child{c,1};
                            child.parent = [wf.frame, wf.index];
                            wavefronts{1,k}{child.index,1} = child;
                        end
                        wavefronts        = backfill_wf_path(wavefronts, wf, struct());
                        wavefront_dynamics = update_wf_database(wavefront_dynamics, wf, wavefronts, threshold);
                        wf_terminated = true;
                        break;

                    case 'MERGED'
                        wf.path           = path_buf(1:path_len, :);
                        wf.length         = len_buf(1:path_len);
                        wf.cause_of_death = 'MERGED';
                        wf.child = [wf_child.frame, wf_child.index];
                        wf.mate  = [wf_mate.frame,  wf_mate.index];

                        % Register child's parents
                        wf_child.parent = [wf.frame, wf.index; wf_mate.frame, wf_mate.index];
                        wavefronts{1,wf_child.frame}{wf_child.index,1} = wf_child;

                        % Register mate's bookkeeping
                        wf_mate.mate           = [wf.frame, wf.index];
                        wf_mate.child          = wf.child;
                        wf_mate.cause_of_death = 'MERGED';
                        wavefronts{1,wf_mate.frame}{wf_mate.index,1} = wf_mate;

                        extra.child = wf.child;
                        extra.mate  = wf.mate;
                        wavefronts        = backfill_wf_path(wavefronts, wf, extra);
                        wavefront_dynamics = update_wf_database(wavefront_dynamics, wf, wavefronts, threshold);
                        wf_terminated = true;
                        break;
                end

            else
                % Empty frame — wavefront expired due to gap
                wf.path           = path_buf(1:path_len, :);
                wf.length         = len_buf(1:path_len);
                wf.cause_of_death = 'EXPIRED';
                wavefronts        = backfill_wf_path(wavefronts, wf, struct());
                wavefront_dynamics = update_wf_database(wavefront_dynamics, wf, wavefronts, threshold);
                wf_terminated = true;
                break;
            end
        end

        % Wavefront alive through the last frame — record it
        if ~wf_terminated
            wf.path           = path_buf(1:path_len, :);
            wf.length         = len_buf(1:path_len);
            wf.cause_of_death = 'EXPIRED';
            wavefronts        = backfill_wf_path(wavefronts, wf, struct());
            wavefront_dynamics = update_wf_database(wavefront_dynamics, wf, wavefronts, threshold);
        end
    end
end

catch ME
    disp(ME);

end

end % end of main function


%% Local helper — writes terminal wf state back to every frame in wf.path
function wavefronts = backfill_wf_path(wavefronts, wf, extra)
if isempty(wf.path)
    return;
end
extra_fields = fieldnames(extra);
num_links = size(wf.path, 1);
for p = 1:num_links
    frame = wf.path(p, 1);
    index = wf.path(p, 2);
    entry = wavefronts{1,frame}{index,1};
    entry.cause_of_death = wf.cause_of_death;
    entry.lifespan       = wf.lifespan;
    entry.path           = wf.path;
    entry.length         = wf.length;
    entry.birthday       = wf.birthday;
    entry.neighbors      = wf.neighbors;
    for f = 1:numel(extra_fields)
        entry.(extra_fields{f}) = extra.(extra_fields{f});
    end
    wavefronts{1,frame}{index,1} = entry;
end
end





