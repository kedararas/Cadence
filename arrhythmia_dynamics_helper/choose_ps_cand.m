function [ps_status, ps_cand] = choose_ps_cand(ps, current_ps_list, next_ps_list, threshold, expected_pos)
% choose_ps_cand  Find the best matching PS candidate in the next frame.
%
%   ps              : current Phase_Singularity object being tracked
%   current_ps_list : cell array of PS in the current frame
%   next_ps_list    : cell array of PS candidates in the next frame
%   threshold       : max spatial distance (pixels) to accept a match
%   expected_pos    : [x, y] predicted position from velocity extrapolation
%                     (optional — uses ps.location if omitted)
%
% Matching enforces two constraints:
%   1. Spatial proximity — candidate must be within threshold of expected_pos
%   2. Charge conservation — candidate must have the same topological charge

% Early exit if either frame has no PS
if isempty(current_ps_list) || isempty(next_ps_list)
    ps_status = 'EXPIRED';
    ps_cand   = [];
    return;
end

% Use velocity-extrapolated position if supplied, otherwise use current location
if nargin >= 5 && ~isempty(expected_pos)
    query_pos = expected_pos;
else
    first_list = cell2mat(cellfun(@(p) p.location, current_ps_list, 'UniformOutput', false));
    query_pos  = first_list(ps.index, :);
end

% Extract next-frame locations — reshape to n×2 regardless of next_ps_list orientation
second_list = cell2mat(cellfun(@(p) p.location, next_ps_list(:), 'UniformOutput', false));

% Compute distances from query position to all candidates — 1×m only
dist_row = pdist2(query_pos, second_list);

% Charge conservation constraint — mask out candidates with wrong charge.
% next_ps_list(:).' forces 1×n so charge_mask always matches dist_row shape.
charge_mask = cellfun(@(p) p.charge == ps.charge, next_ps_list(:).');
dist_row(~charge_mask) = Inf;

% Find closest charge-compatible candidate
[min_dist, min_idx] = min(dist_row);

if min_dist < threshold
    ps_status = 'ALIVE';
    ps_cand   = next_ps_list{min_idx};
else
    ps_status = 'EXPIRED';
    ps_cand   = [];
end
