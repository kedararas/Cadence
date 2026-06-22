function relative_dist = get_frechet_distance_single(first_list, second_list)
% Compute discrete Fréchet distance between corresponding frames of two
% wavefront sequences. Lists may have different lengths — the shorter one's
% last frame is reused for the remaining entries (hold-last semantics).

relative_dist = [];
if isempty(first_list) || isempty(second_list)
    return;
end

num_first  = size(first_list,  1);
num_second = size(second_list, 1);
num_wf     = max(num_first, num_second);

relative_dist = nan(num_wf, 1);
first_wf  = [];
second_wf = [];

for i = 1:num_wf
    if i <= num_first
        first_wf  = first_list{i};
    end
    if i <= num_second
        second_wf = second_list{i};
    end
    if ~isempty(first_wf) && ~isempty(second_wf)
        relative_dist(i) = DiscreteFrechetDist(first_wf, second_wf);
    end
end

end
