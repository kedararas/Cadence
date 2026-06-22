function [repeat_wf, f_dist, e_dist] = check_for_wf_repeatability(wf_dynamics, wf, wavefronts, threshold)
% check_for_wf_repeatability  Test whether wf recurs as an earlier tracked wavefront.
%
% Returns the index into wf_dynamics.wf_multiplicity of the matching wavefront
% (repeat_wf > 0), or 0 if no match found. f_dist and e_dist are the median
% Fréchet and Euclidean distances of the best match, or 0 if no match.

repeat_wf = 0;
f_dist    = 0;
e_dist    = 0;

wf_list = wf_dynamics.wf_multiplicity;
if isempty(wf_list) || wf.lifespan(1,2) < threshold(1,8) || round(median(wf.length)) < 6*threshold(1,7)
    return;
end

% wf.length is a column vector (one entry per path frame) set during tracking.
% Use it directly instead of re-reading location matrix sizes from wavefronts.
path_length      = size(wf.path, 1);
first_list_length = wf.length;  % already computed, no cell lookup needed

% Build the reference location sequence for this wavefront — done once per call.
first_list = cell(path_length, 1);
for i = 1:path_length
    first_list{i} = wavefronts{1, wf.path(i,1)}{wf.path(i,2), 1}.location;
end

num_candidates = size(wf_list, 1);
for j = 1:num_candidates

    % --- Pre-filter 1: wavefront size similarity ---
    % Use the stored wf.length vector — avoids traversing location matrices.
    wf_to_compare   = wavefronts{1, wf_list(j,2)}{wf_list(j,3), 1};
    second_list_length = wf_to_compare.length;  % already computed during tracking

    len_diff = abs(mean(first_list_length) - mean(second_list_length));
    if len_diff > threshold(1,1) * mean(second_list_length)
        continue;
    end

    % --- Pre-filter 2: lifespan similarity ---
    time = abs(wf.lifespan(1,2) - wf_list(j,5));
    if time > threshold(1,2) * wf_list(j,5)
        continue;
    end

    % --- Both filters passed: build candidate location sequence ---
    wf_path_length = size(wf_to_compare.path, 1);
    second_list    = cell(wf_path_length, 1);
    for k = 1:wf_path_length
        second_list{k} = wavefronts{1, wf_to_compare.path(k,1)}{wf_to_compare.path(k,2), 1}.location;
    end

    % --- Shape comparison ---
    % Fréchet distance: captures geometric similarity of the wavefront trajectory
    frech_dist = get_frechet_distance_single(first_list, second_list);

    % Euclidean distance: point-to-point median distance between corresponding frames.
    % Computed directly — avoids the O(N²) DTW dynamic programming that get_dtw_distance_single
    % also computes but whose output (dtw_dist) was never used by the acceptance criterion.
    num_wf     = max(path_length, wf_path_length);
    eucl_dist  = nan(num_wf, 1);
    for q = 1:num_wf
        fw = first_list{min(q, path_length)};
        sw = second_list{min(q, wf_path_length)};
        if ~isempty(fw) && ~isempty(sw)
            eucl_dist(q) = median(diag(pdist2(fw, sw, 'euclidean')));
        end
    end

    f_dist = median(frech_dist);
    e_dist = median(eucl_dist);

    if f_dist <= 2*threshold(1,7) && e_dist <= threshold(1,7)
        repeat_wf = j;
        return;
    end
end

end
