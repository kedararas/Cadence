function wavefronts = create_wavefront_objects(wavefront_list, threshold)
%CREATE_WAVEFRONT_OBJECTS  Wrap per-frame segment lists as Wavefront objects.
%
%   wavefronts = create_wavefront_objects(wavefront_list, threshold)
%
%   Converts the raw per-frame segmented wavefront coordinates into arrays of
%   Wavefront handle objects, dropping any segment shorter than threshold
%   pixels (too small to track reliably). Each kept object is stamped with its
%   frame index and an intra-frame counter and carries its location matrix and
%   length. Consumed by extract_wavefront_dynamics for tracking.
%
%   Inputs
%     wavefront_list  1 x num_frames cell; each cell is itself a cell array of
%                     2 x L coordinate matrices (one per detected segment).
%     threshold       minimum segment length (pixels) to keep.
%
%   Output
%     wavefronts      1 x num_frames cell; each non-empty cell holds an N x 1
%                     cell of Wavefront objects for that frame.
wavefronts = {};
if isempty(wavefront_list)
    return;
end

num_frames = numel(wavefront_list);
wavefronts = cell(1, num_frames);

for i = 1:num_frames
    wf_list = wavefront_list{i};
    if isempty(wf_list)
        continue;
    end

    num_wf      = numel(wf_list);
    object_list = cell(num_wf, 1);  % pre-allocate; trim unused slots after loop
    wf_counter  = 0;

    for j = 1:num_wf
        location  = wf_list{j}';
        wf_length = size(location, 1);
        if wf_length >= threshold
            wf_counter = wf_counter + 1;
            wf          = Wavefront(i, wf_counter);
            wf.location = location;
            wf.length   = wf_length;
            object_list{wf_counter} = wf;
        end
    end

    wavefronts{i} = object_list(1:wf_counter, 1);  % ,1 preserves N×1 column orientation
    % Indexing with a row vector (1:N) would return a 1×N row cell; adding
    % the explicit column subscript keeps the result as an N×1 column cell.
    % A row cell propagates into cellfun/accumarray and causes a subscript
    % dimension mismatch ('MATLAB:accumarray:badSizeValInputMatInd').
end


