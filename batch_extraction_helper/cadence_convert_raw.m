function d = cadence_convert_raw(entry, fs_default, logf)
%CADENCE_CONVERT_RAW  Read one raw recording into a cmos_all_data struct.
%
%   cmos_all_data = cadence_convert_raw(entry)
%   cmos_all_data = cadence_convert_raw(entry, fs_default, logf)
%
%   Headless port of the per-recording conversion in
%   Cadence_Data_Conversion.mlapp (CONVERT TO MATLAB DATA FORMAT).  `entry`
%   is one element of cadence_raw_plan(folder).  Each camera file becomes
%   CAM<n> (double, rows x cols x frames) with a CAM<n>_image thumbnail:
%     .tif   tiffreadVolume; thumbnail = first frame rescaled to [0 255]
%     .gsh   GSDconverter_BVW (SciMedia .gsh/.gsd); camera 1 also supplies
%            acqFreq, and analog1/analogTime when a regular pacing channel is
%            found in the analog block.
%   num_files is the camera count.  acqFreq falls back to fs_default (1000 Hz
%   when not given -- the app's Sampling Frequency field default) if the raw
%   files do not carry it (.tif always, .gsh never).
%
%   The struct is NOT schema-checked here; the caller runs
%   schema_gate(d, "raw", name) and decides what to do with a failure, as the
%   Signal Conditioning app does before conditioning.
%
%   See also cadence_raw_plan, cadence_condition_data, cadence_batch_pipeline.

    if nargin < 2 || isempty(fs_default) || ~isfinite(fs_default) || fs_default <= 0
        fs_default = 1000;
    end
    if nargin < 3 || isempty(logf), logf = @(s) fprintf('%s\n', s); end

    cams = entry.cam_names;
    ncam = numel(cams);
    d = struct();

    if strcmp(entry.ext, 'tif')
        for c = 1:ncam
            cam_field = sprintf('CAM%d', c);
            img_field = sprintf('CAM%d_image', c);
            vol = double(tiffreadVolume(fullfile(entry.folder, cams{c})));
            d.(cam_field) = vol;

            first_frame = vol(:,:,1);
            lo = min(first_frame(:));
            hi = max(first_frame(:));
            if hi > lo
                d.(img_field) = round(rescale_sat(first_frame, [lo hi], [0 255]));
            else
                % Flat/blank first frame -- avoid divide-by-zero in rescale_sat
                d.(img_field) = zeros(size(first_frame));
            end
            logf(sprintf('%s <- %s (%dx%dx%d)', cam_field, cams{c}, size(vol,1), size(vol,2), size(vol,3)));
        end
    else   % gsh / gsd
        for c = 1:ncam
            d = GSDconverter_BVW(d, entry.folder, cams{c}, c);
            cam_field = sprintf('CAM%d', c);
            if isfield(d, cam_field)
                logf(sprintf('%s <- %s (%dx%dx%d)', cam_field, cams{c}, ...
                    size(d.(cam_field),1), size(d.(cam_field),2), size(d.(cam_field),3)));
            else
                logf(sprintf('%s <- %s: converter returned no data', cam_field, cams{c}));
            end
        end
    end
    d.num_files = ncam;

    if ~isfield(d, 'acqFreq')
        d.acqFreq = fs_default;
        logf(sprintf('acqFreq not in raw files; using %g Hz', fs_default));
    end
end
