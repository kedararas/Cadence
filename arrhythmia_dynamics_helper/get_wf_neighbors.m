function wf_neighbors = get_wf_neighbors(wf, wf_list, neighbor_threshold, debug_mode)
%GET_WF_NEIGHBORS  Find wavefronts spatially adjacent to a given wavefront.
%
%   wf_neighbors = get_wf_neighbors(wf, wf_list, neighbor_threshold, debug_mode)
%
%   Returns the subset of wf_list whose pixel locations come within
%   neighbor_threshold of any pixel of wf. Distances are computed in one
%   vectorised pdist2 call over all candidate locations (tagged back to their
%   source wavefront via accumarray), so the cost is independent of how the
%   candidates are split across objects. Used during wavefront tracking to
%   restrict candidate matching to a local neighbourhood.
%
%   Inputs
%     wf                  Wavefront handle object (uses wf.location, M x 2).
%     wf_list             cell array of candidate Wavefront objects.
%     neighbor_threshold  max pixel distance to count as a neighbour.
%     debug_mode          flag for optional diagnostic plotting/printing.
%
%   Output
%     wf_neighbors        cell array (subset of wf_list) within threshold;
%                         [] if wf_list is empty or wf has no location.
wf_neighbors = [];

% isempty(wf) never fires for a handle object — guard on wf.location instead.
if isempty(wf_list) || isempty(wf.location)
    return;
end

% Force column-cell orientation.  create_wavefront_objects returns N×1 cells,
% but defensive coercion here prevents accumarray subscript-dimension errors
% if any caller passes a row cell.
% Background: cellfun inherits the cell's orientation, so a 1×N row cell
% produces a 1×N row wf_ids.  accumarray interprets a 1×N row subs as ONE
% multi-dimensional subscript (1 row, N dims) and expects val to have exactly
% 1 element → 'MATLAB:accumarray:badSizeValInputMatInd'.
wf_list = wf_list(:);   % N×1 column cell, always
num_wf  = numel(wf_list);

% Vectorized: concatenate all candidate locations into one matrix,
% track which wavefront each row belongs to, then call pdist2 once.
wf_sizes  = cellfun(@(c) size(c.location, 1), wf_list);  % N×1
wf_ids    = repelem((1:num_wf)', wf_sizes(:));
wf_ids    = wf_ids(:);   % force N×1 column — when num_wf==1, (1:1)' is treated
                         % as a scalar by repelem and returns a 1×n row vector

% vertcat is orientation-agnostic and always stacks location matrices
% vertically, unlike cell2mat(cellfun(...)) which would concatenate
% horizontally if the intermediate cell were a row cell.
locs_cell = cellfun(@(c) c.location, wf_list, 'UniformOutput', false);
all_locs  = vertcat(locs_cell{:});                        % sum(wf_sizes) × 2

% One pdist2 call: rows = wf query pts, cols = all candidate pts
D = pdist2(wf.location, all_locs);                  % n_query × n_all
min_per_col = min(D, [], 1);                        % min over query pts → 1 × n_all

% Per-wavefront minimum distance — both subs and val must be column vectors.
per_wf_min = accumarray(wf_ids, min_per_col(:), [num_wf, 1], @min);

% Build [min_dist, wf_index, wf_size] table; filter before sorting
euclidean_relative_distance = [per_wf_min, (1:num_wf)', wf_sizes(:)];
euclidean_relative_distance(euclidean_relative_distance(:,1) > neighbor_threshold, :) = [];

if ~isempty(euclidean_relative_distance)
    wf_neighbors = sortrows(euclidean_relative_distance, 1);
end

if debug_mode && ~isempty(wf_list)
    % Cache colors outside loop; reuse figure if one is already open
    col_black = [0 0 0];
    col_blue  = [0 0 1];
    col_red   = [1 0 0];

    figure; hold on; axis fill; axis off;

    for wf_count = 1:num_wf
        loc = wf_list{wf_count}.location;
        if ~isempty(loc)
            plot(loc(:,1), loc(:,2), 'Color', col_black, 'LineWidth', 3);
        end
    end

    plot(wf.location(:,1), wf.location(:,2), 'Color', col_blue, 'LineWidth', 3);

    for i = 1:size(wf_neighbors, 1)
        loc = wf_list{wf_neighbors(i,2)}.location;
        if ~isempty(loc)
            plot(loc(:,1), loc(:,2), 'Color', col_red, 'LineWidth', 3);
        end
    end
end

