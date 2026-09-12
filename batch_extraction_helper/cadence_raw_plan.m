function plan = cadence_raw_plan(folder)
%CADENCE_RAW_PLAN  List the recordings held in one folder of raw camera files.
%
%   plan = cadence_raw_plan(folder)
%
%   Headless port of the folder-planning step of Cadence_Data_Conversion.mlapp
%   (CONVERT TO MATLAB DATA FORMAT).  A folder normally holds ONE recording
%   whose files are its 1-4 cameras (.tif volumes, or SciMedia .gsh/.gsd
%   pairs).  For backward compatibility a folder may hold MULTIPLE
%   recordings: when EVERY file of an extension ends in a camera letter A-D
%   immediately before the extension (003-RC-trans-60bpm0311-004A.gsh), the
%   letter is the camera and the rest of the name is the recording.  Files
%   that do not follow that convention are the cameras of a single recording
%   named after the folder, in directory (alphabetical) order.
%
%   Output: struct array, empty when the folder has no .tif / .gsh files.
%     .folder      the folder
%     .mat_base    base name of the .mat to write (<folder name> for a
%                  single recording, the parsed base otherwise)
%     .ext         'tif' | 'gsh'
%     .cam_names   cellstr of file names in CAM1..CAMn order
%
%   See also cadence_convert_raw, cadence_batch_pipeline.

    folder = char(folder);
    folder = regexprep(folder, '[\\/]+$', '');
    parts = strsplit(folder, {'/', '\'});
    folder_name = parts{end};

    plan = struct('folder', {}, 'mat_base', {}, 'ext', {}, 'cam_names', {});
    exts = {'tif', 'gsh'};
    for e = 1:numel(exts)
        ext = exts{e};
        listing = dir(fullfile(folder, ['*.' ext]));
        listing = listing(~[listing.isdir]);
        listing = listing(~startsWith({listing.name}, '.'));   % ._ resource forks on ExFAT/SMB
        if isempty(listing), continue; end
        fnames = {listing.name};

        bases = cell(size(fnames));
        for k = 1:numel(fnames)
            [~, bases{k}, ~] = fileparts(fnames{k});
        end

        % Camera-letter convention present only if EVERY file ends in A-D.
        hasLetter = ~cellfun('isempty', regexp(bases, '[A-Da-d]$', 'once'));
        if all(hasLetter)
            grp    = cellfun(@(s) s(1:end-1), bases, 'UniformOutput', false);
            camIdx = cellfun(@(s) double(upper(s(end))) - double('A') + 1, bases);
            [uniqBase, ~, ic] = unique(grp, 'stable');
            for g = 1:numel(uniqBase)
                sel = find(ic == g);
                [~, ord] = sort(camIdx(sel));   % A,B,C,D order
                sel = sel(ord);
                p = numel(plan) + 1;
                plan(p).folder    = folder;
                plan(p).mat_base  = uniqBase{g};
                plan(p).ext       = ext;
                plan(p).cam_names = fnames(sel);   % compact, no gaps
            end
        else
            % No convention: all files are cameras of one recording.
            p = numel(plan) + 1;
            plan(p).folder    = folder;
            plan(p).mat_base  = folder_name;
            plan(p).ext       = ext;
            plan(p).cam_names = fnames;
        end
    end

    % A folder that collapses to a single recording keeps its legacy
    % <folder>.mat name; multiple recordings keep their parsed bases.
    if numel(plan) == 1
        plan(1).mat_base = folder_name;
    end
end
