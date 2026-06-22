function [wf_status, wf_cand, wf_mate, wf_child] = choose_wf_cand(wf, current_wf_list, next_wf_list, next_wf_black_list, neighbor_threshold, debug_mode)

% Cache frequently reused values
frech_thresh = neighbor_threshold * 4;
wf_len       = wf.length(end, 1);

%% Step 1: Compute Frechet matrix once for this frame pair, then look up rows as needed

% Extract locations from both lists once
current_locs = cellfun(@(c) c.location, current_wf_list, 'UniformOutput', false);
next_locs    = cellfun(@(c) c.location, next_wf_list,    'UniformOutput', false);
frech_matrix = get_frechet_distance(current_locs, next_locs); % n×m, computed once

wf_future_neighbors = get_wf_neighbors(wf, next_wf_list, neighbor_threshold, debug_mode);
wf_similars         = get_similar_wf(wf.index, next_wf_list, frech_matrix, frech_thresh);

% Convenience flag — avoids repeating the same multi-part empty check below.
% Note: all neighbor table accesses use && (short-circuit) so the right-hand
% side is only evaluated when the guard is true.
has_wf_neighbors = ~isempty(wf_future_neighbors);

% Evaluate neighbor guard once — used by both closest and second_closest blocks
has_current_neighbor = ~isempty(wf.neighbors) && ...
    ~isempty(wf.neighbors{end,1}) && ...
    wf.neighbors{end,1}(1,3) == wf.frame;

if has_current_neighbor
    closest_wf                  = current_wf_list{wf.neighbors{end,1}(1,2), 1};
    closest_wf_future_neighbors = get_wf_neighbors(closest_wf, next_wf_list, neighbor_threshold, debug_mode);
    closest_wf_similars         = get_similar_wf(closest_wf.index, next_wf_list, frech_matrix, frech_thresh);
    has_closest_neighbors       = ~isempty(closest_wf_future_neighbors);
else
    closest_wf                  = [];
    closest_wf_future_neighbors = [];
    closest_wf_similars         = [];
    has_closest_neighbors       = false;
end

%% Scenario 1: wf and its closest neighbor merged into one future wavefront

merged = 0;
second = 0;

if has_current_neighbor && has_closest_neighbors && has_wf_neighbors
    if closest_wf_future_neighbors(1,2) == wf_future_neighbors(1,2) ...
            && closest_wf_future_neighbors(1,1) < 1 && wf_future_neighbors(1,1) < 1
        num_wf_neigh         = sum(wf_future_neighbors(:,1) < 1);
        num_closest_wf_neigh = sum(closest_wf_future_neighbors(:,1) < 1);
        if num_wf_neigh == 1 && num_closest_wf_neigh == 1
            merged = 1;
        elseif num_wf_neigh > 1 || num_closest_wf_neigh > 1
            if abs((wf_len + closest_wf.length(end,1)) - wf_future_neighbors(1,3)) < 0.2*wf_future_neighbors(1,3)
                merged = 1;
            end
        end
        if ~isempty(closest_wf_similars) && ~isempty(wf_similars) ...
                && closest_wf_similars(1,2) ~= wf_similars(1,2)
            merged = 0;
        end
    end
end

% Only compute second_closest if closest-based merge check failed
if ~merged && has_current_neighbor && size(wf.neighbors{end,1}, 1) > 1
    second_closest_wf                  = current_wf_list{wf.neighbors{end,1}(2,2), 1};
    second_closest_wf_future_neighbors = get_wf_neighbors(second_closest_wf, next_wf_list, neighbor_threshold, debug_mode);
    second_closest_wf_similars         = get_similar_wf(second_closest_wf.index, next_wf_list, frech_matrix, frech_thresh);
    has_second_closest_neighbors       = ~isempty(second_closest_wf_future_neighbors);

    if has_second_closest_neighbors && has_wf_neighbors
        if second_closest_wf_future_neighbors(1,2) == wf_future_neighbors(1,2) ...
                && second_closest_wf_future_neighbors(1,1) < 1 && wf_future_neighbors(1,1) < 1
            num_wf_neigh                = sum(wf_future_neighbors(:,1) < 1);
            num_second_closest_wf_neigh = sum(second_closest_wf_future_neighbors(:,1) < 1);
            if num_wf_neigh == 1 && num_second_closest_wf_neigh == 1
                merged = 1; second = 1;
            elseif num_wf_neigh > 1 || num_second_closest_wf_neigh > 1
                if abs((wf_len + second_closest_wf.length(end,1)) - wf_future_neighbors(1,3)) < 0.2*wf_future_neighbors(1,3)
                    merged = 1; second = 1;
                end
            end
            if ~isempty(second_closest_wf_similars) && ~isempty(wf_similars) ...
                    && second_closest_wf_similars(1,2) ~= wf_similars(1,2)
                merged = 0;
            end
        end
    end
else
    second_closest_wf = [];
end

% merged=1 is only reachable through blocks that already verified
% has_wf_neighbors, so wf_future_neighbors(1,2) is safe here.
if merged
    wf_status = 'MERGED';
    wf_child  = next_wf_list{wf_future_neighbors(1,2), 1};
    wf_mate   = second_closest_wf;
    if ~second
        wf_mate = closest_wf;
    end
    wf_cand = [];
    return;
end

%% Scenario 2: wavefront fragmented into two child wavefronts

fragmented = 0;

% Guard wf_future_neighbors before indexing — using & on an empty matrix
% throws 'MATLAB:badsubscript' because both sides of & are always evaluated.
if has_wf_neighbors
    num_close_neighbors = sum(wf_future_neighbors(:,1) < 1);
else
    num_close_neighbors = 0;
end

if num_close_neighbors >= 2
    % num_close_neighbors >= 2 guarantees at least 2 rows in wf_future_neighbors
    if (~isempty(wf_similars)         && wf_similars(1,1)         < neighbor_threshold) || ...
       (~isempty(closest_wf_similars) && closest_wf_similars(1,1) < neighbor_threshold)
        fragmented = 0; % similar to itself — not a true fragmentation
    else
        wf_n1       = next_wf_list{wf_future_neighbors(1,2), 1}.length(1,1);
        wf_n2       = next_wf_list{wf_future_neighbors(2,2), 1}.length(1,1);
        frag_length = wf_n1 + wf_n2;
        if frag_length <= wf_len || (wf_len > wf_n1 && wf_len > wf_n2)
            fragmented = 1;
        end
    end
end

if fragmented
    wf_status     = 'FRAGMENTED';
    wf_child{1,1} = next_wf_list{wf_future_neighbors(1,2), 1};
    wf_child{2,1} = next_wf_list{wf_future_neighbors(2,2), 1};
    wf_mate = [];
    wf_cand = [];
    return;
end

%% Scenario 3: wavefront maintained its unique identity

% All three match flags use && (short-circuit) so wf_future_neighbors is
% only indexed when has_wf_neighbors is true.
good_match = ~isempty(wf_similars) && wf_similars(1,1) < neighbor_threshold;
fair_match = ~isempty(wf_similars) && has_wf_neighbors && ...
    wf_similars(1,2) == wf_future_neighbors(1,2) && wf_future_neighbors(1,1) < 1;
size_match = has_wf_neighbors && wf_future_neighbors(1,1) < 1 && ...
    abs(wf_len - wf_future_neighbors(1,3)) < 0.2*wf_len;

if good_match || fair_match || size_match
    if ~isempty(wf_similars)
        wf_cand = next_wf_list{wf_similars(1,2), 1};
        if size(wf_cand.length,1) == 1 && isempty(wf_cand.cause_of_death) && ...
                abs(wf_len - wf_similars(1,3)) < 0.2*wf_len
            wf_status = 'ALIVE';
            wf_mate   = [];
            wf_child  = [];
            return;
        end
    end

    % size([], 1) == 0 so this loop is safely skipped when wf_future_neighbors is empty
    for i = 1:size(wf_future_neighbors, 1)
        wf_cand = next_wf_list{wf_future_neighbors(i,2), 1};
        if size(wf_cand.length, 1) == 1
            wf_status = 'ALIVE';
            wf_mate   = [];
            wf_child  = [];
            return;
        end
    end
end

%% Default: wavefront expired
wf_status = 'EXPIRED';
wf_cand   = [];
wf_mate   = [];
wf_child  = [];
