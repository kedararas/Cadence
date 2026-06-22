function [ps_objects, valid_ps] = create_ps_objects(ps_matrix)
% create_ps_objects  Convert PS point matrices into Phase_Singularity handle objects.
%
% ps_matrix must be a 1×num_frames cell — pass cmos_all_data.ps(1,:), NOT the
% full 3-row ps cell. Each ps_matrix{i} is a numeric matrix with rows:
%   [x, y, charge, wf_dist, wf_index]  (5 columns, from tag_ps_to_wf output)
%
% ps_objects{i}: cell column of Phase_Singularity objects for frame i
% valid_ps{i}:   zero column vector, filled in by update_ps_database as PS are confirmed

ps_objects = {};
valid_ps   = {};

if isempty(ps_matrix)
    return;
end

num_frames = numel(ps_matrix);
ps_objects = cell(1, num_frames);
valid_ps   = cell(1, num_frames);

for i = 1:num_frames
    ps_list = ps_matrix{i};
    if isempty(ps_list)
        continue;
    end

    num_ps      = size(ps_list, 1);
    object_list = cell(num_ps, 1);  % pre-allocate

    for j = 1:num_ps
        row = ps_list(j,:);
        ps            = Phase_Singularity(i, j);
        ps.location   = [row(1), row(2)];   % [x, y]
        ps.charge     = row(3);
        ps.wavefront  = [i, row(5)];        % [frame, wf_index] — col 5 after tag_ps_to_wf refactor
        object_list{j} = ps;
    end

    ps_objects{i} = object_list;
    valid_ps{i}   = zeros(num_ps, 1);
end


