function [ps, ps_count] = count_ps(ps_data, wavefronts, pixel_size_mm)
% count_ps  Detect, deduplicate, and annotate phase singularities per frame.
%
% Wavefront proximity is treated as METADATA, not a hard filter.
% ps{1,i} contains all cleaned PS with wavefront association columns appended.
% ps{2,i} contains all cleaned PS (no wavefront info).
% ps{3,i} contains all raw PS before deduplication.
%
% Downstream code can filter by col 5 > 0 to isolate WF-associated PS,
% or use all rows and treat col 5 as descriptive metadata.
%
% Optional argument pixel_size_mm (mm/pixel) enables physically-grounded
% thresholds. If omitted, thresholds are derived from frame geometry.

num_frames = size(ps_data, 3);
num_rows   = size(ps_data, 1);
num_cols   = size(ps_data, 2);
frame_diag = sqrt(num_rows^2 + num_cols^2);

if nargin >= 3 && ~isempty(pixel_size_mm) && pixel_size_mm > 0
    % Physical thresholds — grounded in tissue scale
    % PS detection clusters span ~1-2 mm; tip association within ~0.5 mm
    dup_threshold      = round(2.0 / pixel_size_mm);  % merge clusters within 2 mm
    wf_assoc_threshold = round(0.5 / pixel_size_mm);  % PS within 0.5 mm of wavefront tip
    fprintf('count_ps: pixel_size=%.3f mm -> dup=%d px, wf_assoc=%d px\n', ...
        pixel_size_mm, dup_threshold, wf_assoc_threshold);
else
    % Geometry-derived fallbacks when pixel size is unknown
    % Deduplication: ~3% of frame diagonal — merges local detection clusters
    % WF association: ~1% of frame diagonal — PS must be close to wavefront tip
    dup_threshold      = max(5,  round(0.03 * frame_diag));
    wf_assoc_threshold = max(3,  round(0.01 * frame_diag));
    fprintf('count_ps: no pixel_size supplied -> dup=%d px, wf_assoc=%d px\n', ...
        dup_threshold, wf_assoc_threshold);
end

% Pre-allocate outputs — 3 rows each
ps       = cell(3, num_frames);
ps_count = zeros(3, num_frames);

for i = 1:num_frames
    ps_frame = round(squeeze(ps_data(:,:,i)));

    % Detect positive and negative topological charges
    [p_rows, p_cols] = find(ps_frame ==  1);
    [n_rows, n_cols] = find(ps_frame == -1);
    ps_points = [p_cols, p_rows,  ones(numel(p_rows),1); ...
                 n_cols, n_rows, -ones(numel(n_rows),1)];

    if isempty(ps_points)
        ps{1,i} = [];
        ps{2,i} = [];
        ps{3,i} = [];
        continue;
    end

    ps{3,i} = ps_points;

    %% Deduplicate: remove PS within dup_threshold pixels of a higher-priority PS
    % Vectorized inner loop: each kept point suppresses all its neighbors at once
    num_ps    = size(ps_points, 1);
    distances = pdist2(ps_points(:,1:2), ps_points(:,1:2));
    keep      = true(num_ps, 1);
    for x = 1:num_ps
        if keep(x)
            neighbors    = distances(x,:) < dup_threshold;
            neighbors(x) = false;           % exclude self
            keep(neighbors) = false;        % suppress all neighbors in one write
        end
    end
    cleaned_ps_points = ps_points(keep, :);

    %% Tag each PS to its nearest wavefront (metadata only — not a filter)
    wf_tagged = tag_ps_to_wf(cleaned_ps_points, wavefronts{2,i}, wf_assoc_threshold);

    ps{1,i} = wf_tagged;
    ps{2,i} = cleaned_ps_points;

    ps_count(1,i) = sum(wf_tagged(:,5) > 0);       % PS with a wavefront association
    ps_count(2,i) = size(cleaned_ps_points, 1);     % all cleaned PS
    ps_count(3,i) = size(ps_points, 1);             % raw PS before dedup
end
