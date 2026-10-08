function [stamp, raw] = cadence_acquisition_time(entry)
%CADENCE_ACQUISITION_TIME  The camera's own acquisition time stamp of a recording.
%
%   stamp = cadence_acquisition_time(plan_entry)        % one entry of cadence_raw_plan
%   stamp = cadence_acquisition_time('/path/CAM1.gsh')  % or a header / .tif file
%   [stamp, raw] = ...                                   % raw = the text as found
%
%   Reads ONLY the header, never the image data, so it takes milliseconds per
%   recording and needs no conversion:
%     .gsh (BV Workbench)  the line  "Date created : 12/8/2025 8:45:39 PM"
%     .gsh (MiCAM)         the line  "AcquisitionDate=..."
%     .tif                 the TIFF DateTime tag of the first file
%   stamp is "yyyy-mm-dd HH:MM:SS" (a string), or "" when the header has no
%   time stamp or it cannot be parsed (raw then still holds the text found).
%   This is the acquisition computer's clock; the Acquired column of the batch
%   tables is the nominal time written in the file label.

    stamp = "";  raw = "";
    if isstruct(entry)
        if isempty(entry.cam_names), return; end
        file = fullfile(entry.folder, entry.cam_names{1});
    else
        file = char(entry);
    end
    if ~isfile(file), return; end
    [~, ~, ext] = fileparts(file);
    try
        if strcmpi(ext, '.gsh')
            fid = fopen(file, 'r');
            if fid < 0, return; end
            txt = fread(fid, 4096, '*char')';      % the header is a few hundred bytes of text
            fclose(fid);
            t = regexp(txt, 'Date created\s*:\s*([^\r\n]+)', 'tokens', 'once');
            if isempty(t), t = regexp(txt, 'AcquisitionDate\s*=\s*([^\r\n]+)', 'tokens', 'once'); end
            if isempty(t), return; end
            raw = string(strtrim(t{1}));
        else
            info = imfinfo(file);
            if ~isfield(info, 'DateTime') || isempty(info(1).DateTime), return; end
            raw = string(strtrim(info(1).DateTime));
        end
    catch
        return;
    end
    fmts = {'M/d/yyyy h:mm:ss a', 'yyyy/MM/dd HH:mm:ss', 'yyyy-MM-dd HH:mm:ss', 'yyyy:MM:dd HH:mm:ss', ...
            'MM/dd/yyyy HH:mm:ss', 'yyyy/MM/dd HH:mm', 'd-MMM-yyyy HH:mm:ss'};
    for k = 1:numel(fmts)
        try
            d = datetime(raw, 'InputFormat', fmts{k}, 'Locale', 'en_US');
            if ~isnat(d), stamp = string(d, 'yyyy-MM-dd HH:mm:ss'); return; end
        catch
        end
    end
end
