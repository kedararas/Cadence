function wf_cand_list = get_similar_wf(wf_index, next_wf_list, frech_matrix, threshold)
% Returns sorted Frechet-similar candidates from next_wf_list for wavefront wf_index.
%
%   wf_index    : row index into frech_matrix identifying the current wavefront
%   next_wf_list: cell array of candidate wavefront objects in the next frame
%   frech_matrix: pre-computed n×m Frechet distance matrix
%                 (rows = current frame wfs, cols = next frame wfs)
%                 Compute once per frame pair via get_frechet_distance and pass here.
%   threshold   : maximum Frechet distance to accept as a candidate

wf_cand_list = [];

if isempty(next_wf_list) || isempty(frech_matrix)
    return;
end

% Extract and sort the single row for this wavefront
[sorted_dist, sorted_index] = sort(frech_matrix(wf_index, :));

% Filter by threshold
mask = sorted_dist <= threshold;
if ~any(mask)
    return;
end

sorted_dist  = sorted_dist(mask)';
sorted_index = sorted_index(mask)';

% Build output table: [frechet_dist, wf_index_in_next, wf_size]
wf_sizes     = cellfun(@(c) size(c.location, 1), next_wf_list(sorted_index));
wf_cand_list = [sorted_dist, sorted_index, wf_sizes(:)];
