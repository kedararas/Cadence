function relative_dist = get_frechet_distance(first_list, second_list)
%GET_FRECHET_DISTANCE  Pairwise discrete Frechet distances between wavefront sets.
%
%   relative_dist = get_frechet_distance(first_list, second_list)
%
%   Builds the full M x N matrix of discrete Frechet distances between every
%   wavefront curve in first_list and every curve in second_list, by calling
%   DiscreteFrechetDist on each pair. Computed once per frame pair and then
%   indexed by get_similar_wf during candidate matching. Pairs that return no
%   distance are recorded as -Inf so they sort as "most dissimilar".
%
%   Inputs
%     first_list   M x 1 cell, each cell a P x 2 wavefront location matrix.
%     second_list  N x 1 cell, each cell a Q x 2 wavefront location matrix.
%
%   Output
%     relative_dist  M x N double; relative_dist(j,k) is the Frechet distance
%                    between first_list{j} and second_list{k} (NaN-initialised,
%                    -Inf where no distance was returned); [] if either list
%                    is empty.
    relative_dist = [];
    if isempty(first_list) || isempty(second_list)
        return; %nothing to compare
    end
    
    first_size = size(first_list,1);
    second_size = size(second_list,1);
    
    relative_dist = nan(first_size, second_size);
    for j=1:first_size
        first_wf_j = first_list{j,1};
        for k=1:second_size
            second_wf_k = second_list{k,1};
            [frech_dist, ~] = DiscreteFrechetDist(first_wf_j, second_wf_k);
            if ~isempty(frech_dist)
                relative_dist(j,k) = frech_dist;
            else
                relative_dist(j,k) = -Inf;
            end
        end
    end

end

