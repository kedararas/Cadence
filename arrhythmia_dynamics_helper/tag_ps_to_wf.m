function ps_locations = tag_ps_to_wf(ps_points, wavefronts, threshold)
% tag_ps_to_wf  Annotate each PS with its nearest wavefront (metadata only).
%
% Returns ALL ps_points with two extra columns appended:
%   col 4: distance to nearest wavefront endpoint (Inf if no wavefronts present)
%   col 5: index of nearest wavefront (0 if distance > threshold)
%
% No PS are discarded. Callers can filter by col 5 > 0 to isolate
% WF-associated PS, or use all rows and treat col 5 as descriptive metadata.

num_ps = size(ps_points, 1);

% Default: all PS get Inf distance and wavefront index 0 (not associated)
ps_locations = [ps_points, inf(num_ps,1), zeros(num_ps,1)];

if isempty(wavefronts) || isempty(ps_points)
    return;
end

num_wf = numel(wavefronts);

%% Build wf endpoint matrix — start and end point of each wavefront
% Interleaved: [wf1_start; wf1_end; wf2_start; wf2_end; ...]
wf_points = zeros(num_wf*2, 3);
for i = 1:num_wf
    wf = wavefronts{i};
    wf_points(2*i-1, :) = [wf(1,1),   wf(2,1),   i];
    wf_points(2*i,   :) = [wf(1,end),  wf(2,end), i];
end

%% Compute distances and pre-sort rows and columns once
relative_dist = pdist2(wf_points(:,1:2), ps_points(:,1:2));  % (2*num_wf) × num_ps

% Pre-sort rows: for each wf endpoint, PS candidates in distance order
[ps_dist_sorted, ps_idx_sorted] = sort(relative_dist, 2);

% Pre-sort columns: for each PS, wf endpoints in distance order
[~, wf_idx_sorted] = sort(relative_dist, 1);

%% Greedy mutual nearest matching
% tagged: prevents one PS from being claimed by multiple wavefronts
tagged = false(num_ps, 1);

for j = 1:num_wf*2
    for k = 1:num_ps
        cand_ps = ps_idx_sorted(j, k);
        if tagged(cand_ps)
            continue;
        end

        % Check mutual nearest: is this wf endpoint the closest unprocessed
        % endpoint to cand_ps?
        wf_col = wf_idx_sorted(:, cand_ps);
        wf_col(wf_col < j) = [];   % exclude already-processed endpoints

        if ~isempty(wf_col) && wf_col(1) == j
            % Mutual match — annotate with distance and wavefront index
            wf_dist = ps_dist_sorted(j, k);
            wf_idx  = ceil(j / 2);
            ps_locations(cand_ps, 4) = wf_dist;
            ps_locations(cand_ps, 5) = wf_idx * (wf_dist <= threshold);  % 0 if too far
            tagged(cand_ps) = true;
            break;
        end
    end
end
