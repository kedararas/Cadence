function cmos_all_data = GSDconverter_BVW(cmos_all_data, input_path, file_name, camera_counter)
%GSDCONVERTER_BVW  Read a SCIMedia .gsh/.gsd file pair and append to cmos_all_data.
%
%   cmos_all_data = GSDconverter_BVW(cmos_all_data, input_path, file_name, camera_counter)
%
%   Inputs
%     cmos_all_data    – struct being built by the caller
%     input_path       – folder that contains the .gsh / .gsd files
%     file_name        – file name (extension optional)
%     camera_counter   – 1-based camera index (controls which CAM# field is written)
%
%   Outputs
%     cmos_all_data    – struct with CAM<n>, CAM<n>_image, CAM<n>_raw_image appended.
%                        For camera 1, acqFreq, analog1, and analogTime are also added
%                        when a regular pacing signal is detected.
%
%   Supported acquisition software
%     0 = MiCAM  (header uses '=' delimiter, first line contains 'DataName')
%     1 = BV Workbench  (header uses ':' delimiter; requires summary.txt)
%
%   File format notes
%     .gsh  – plain-text header
%     .gsd  – little-endian int16 binary (256-byte reserved + FORM_INFO + AUX_INFO + data)
%     Pixel data starts at byte 972.  Frame 0 is the background image.
%     Analog data follows immediately after all camera frames.

[~, name, ~] = fileparts(file_name);

% ── Read .gsh header ─────────────────────────────────────────────────────────
gsh_path = fullfile(input_path, [name, '.gsh']);
fid_gsh  = fopen(gsh_path, 'r');
if fid_gsh == -1
    error('GSDconverter_BVW:fileNotFound', 'Cannot open header: %s', gsh_path);
end
textInfo    = textscan(fid_gsh, '%s', 'Delimiter', '');
textInfo    = textInfo{1};
fclose(fid_gsh);                              % close only this handle

% Detect acquisition software: MiCAM first line contains 'DataName'
acqVersion = ~contains(textInfo{1}, 'DataName');   % 0 = MiCAM, 1 = BV Workbench

% Defaults (guard against header fields being absent)
acquisitionDate = '';
pageFrames      = 0;
sampleTime      = 1;      % ms
acqFreq         = 1000;   % Hz

if acqVersion == 0
    % ── MiCAM: '=' delimited, arbitrary line order ─────────────────────────
    for k = 1:numel(textInfo)
        ln = textInfo{k};
        eq = strfind(ln, '=');
        if isempty(eq), continue; end
        if     contains(ln, 'AcquisitionDate')
            acquisitionDate = strtrim(ln(eq+1:end));
        elseif contains(ln, 'frame_number')
            pageFrames = str2double(strtrim(ln(eq+1:end)));
        elseif contains(ln, 'sample_time')
            sampleTime = str2double(strtrim(ln(eq+1:end-4)));  % strip trailing 'msec'
            acqFreq    = 1000 / sampleTime;
        end
    end

else
    % ── BV Workbench: ':' delimited, fixed line positions ──────────────────
    if numel(textInfo) >= 2 && contains(textInfo{2}, 'Date created')
        c = strfind(textInfo{2}, ':');
        acquisitionDate = strtrim(textInfo{2}(c(1)+1:end));
    end
    if numel(textInfo) >= 5 && contains(textInfo{5}, 'Number of frames')
        c = strfind(textInfo{5}, ':');
        pageFrames = str2double(strtrim(textInfo{5}(c+1:end)));
    end
    if numel(textInfo) >= 7 && contains(textInfo{7}, 'Frame rate (Hz)')
        c = strfind(textInfo{7}, ':');
        acqFreq    = str2double(strtrim(textInfo{7}(c+1:end)));
        sampleTime = 1000 / acqFreq;
    end
end

% ── BV Workbench: read summary.txt for shutter delay and camera count ─────────
comment = {};
if acqVersion == 1
    txt_path = fullfile(input_path, 'summary.txt');
    fid_txt  = fopen(txt_path, 'r');
    if fid_txt == -1
        warning('GSDconverter_BVW:noSummaryTxt', ...
            'summary.txt not found in %s — camera count unknown.', input_path);
        return;                               % cannot continue without it
    end
    textTxt  = textscan(fid_txt, '%s', 'Delimiter', '');
    textTxt  = textTxt{1};
    fclose(fid_txt);                          % close only this handle

    % Initialize to NaN so we can detect whether the markers were found
    commentStart = NaN;
    commentEnd   = NaN;

    for k = 1:numel(textTxt)
        ln = textTxt{k};
        if contains(ln, 'Active cameras')
            c       = strfind(ln, ':');
            numCams = str2double(strtrim(ln(c+1:end)));
        elseif contains(ln, 'Shutter delay')
            c = strfind(ln, ':');
            % shutterDelay stored for reference; not currently written to struct
        elseif contains(ln, 'Comment')
            commentStart = k;
        elseif contains(ln, 'Tags')
            commentEnd = k;
        end
    end

    % Extract user comment (present only when commentEnd > commentStart + 3)
    if ~isnan(commentStart) && ~isnan(commentEnd) && commentEnd > commentStart + 3
        comment = textTxt(commentStart+2 : commentEnd-2);  %#ok<NASGU>
    end
end

% ── Read .gsd binary ─────────────────────────────────────────────────────────
gsd_path = fullfile(input_path, [name, '.gsd']);
fid_gsd  = fopen(gsd_path, 'r', 'l');         % little-endian
if fid_gsd == -1
    error('GSDconverter_BVW:fileNotFound', 'Cannot open data file: %s', gsd_path);
end

% FORM_INFO (byte 256): pixel geometry
fseek(fid_gsd, 256, 'bof');
numXPixels           = fread(fid_gsd, 1, 'short');
numYPixels           = fread(fid_gsd, 1, 'short');
numXSkippedPixels    = fread(fid_gsd, 1, 'short');
numYSkippedPixels    = fread(fid_gsd, 1, 'short');
numXActualPixels     = fread(fid_gsd, 1, 'short');
numYActualPixels     = fread(fid_gsd, 1, 'short');
fread(fid_gsd, 1, 'short');                   % numFramesGSD (same as pageFrames)

% AUX_INFO (byte 328): analog channel info
fseek(fid_gsd, 328, 'bof');
numAnalogChannels        = fread(fid_gsd, 1, 'short');
analogSamplingMultiplier = fread(fid_gsd, 1, 'short');
fseek(fid_gsd, 338, 'bof');
numAnalogFrames          = fread(fid_gsd, 1, 'short');

% Pixel data starts at byte 972 (background frame + all signal frames)
fseek(fid_gsd, 972, 'bof');
nTotal   = numXPixels * numYPixels * (pageFrames + 1);
rawPx    = fread(fid_gsd, nTotal, 'int16');   % single contiguous read

% Read analog data for camera 1 (follows immediately after pixel block)
rawAnalog = [];
if camera_counter == 1
    nAnalog   = numAnalogChannels * numAnalogFrames * analogSamplingMultiplier;
    rawAnalog = fread(fid_gsd, nAnalog, 'int16');
end

fclose(fid_gsd);                              % close only this handle

% ── Reshape pixel data ────────────────────────────────────────────────────────
% fread fills column-major; x-axis (columns) is fastest, so dim 1 = numXPixels
cmosData = reshape(rawPx, numXPixels, numYPixels, pageFrames + 1);

% Crop ROI (remove skipped border pixels)
cmosData = cmosData(numXSkippedPixels+1 : numXSkippedPixels+numXActualPixels, ...
                    numYSkippedPixels+1  : numYSkippedPixels+numYActualPixels, :);

% Swap x↔y so rows = y-axis, cols = x-axis
cmosData = permute(cmosData, [2, 1, 3]);

% Frame 0 is the background image; frames 1..end are the signal
bgImage  = cmosData(:,:,1);
cmosData = double(cmosData(:,:,2:end));

% Normalize background image to [0, 255].
% Implemented manually to avoid rescale() version incompatibilities and to
% guard against a constant or empty frame (divide-by-zero → all zeros).
bg_min = double(min(bgImage(:)));
bg_max = double(max(bgImage(:)));
if bg_max > bg_min
    bgImageNom = round(255 * (double(bgImage) - bg_min) / (bg_max - bg_min));
else
    bgImageNom = zeros(size(bgImage));
end

% ── Analog pacing signal detection (camera 1 only) ───────────────────────────
analog1    = [];
analogTime = [];

if camera_counter == 1 && ~isempty(rawAnalog)
    % Reshape: rows = (oversampled) time, cols = channels
    analogData = reshape(rawAnalog, numAnalogFrames * analogSamplingMultiplier, []);

    % Downsample to camera frame rate
    analogData = downsample(analogData, analogSamplingMultiplier);
    N_analog   = size(analogData, 1);
    analogTime = (0 : N_analog-1)' / acqFreq;   % seconds

    analog1 = zeros(N_analog, 1);

    for ch_idx = 1:size(analogData, 2)
        ch = double(analogData(:, ch_idx));

        if std(ch) <= 150
            continue;   % too quiet — not a pacing channel
        end

        % detect polarity BEFORE taking abs().
        % Original code applied abs() to all channels first, making isDownward
        % always false and leaving aData undefined when isDownward == 0.
        isDownward = abs(max(ch) - median(ch)) < abs(min(ch) - median(ch));
        if isDownward
            ch = -ch - min(-ch);   % flip so spikes are always upward
        end

        % Threshold at 70 % of spike amplitude above baseline
        delta = 0.7 * (max(ch) - median(ch));
        burst = ch >= (median(ch) + delta);   % logical: 1 during each pacing pulse

        % Find rising edges (start of each burst)
        rising = find(diff([0; burst]) == 1);
        if numel(rising) < 2
            continue;   % not enough beats to assess regularity
        end

        ipi = diff(rising);   % inter-pulse intervals (samples)
        if range(ipi) > 5
            continue;   % too irregular — not a pacing signal
        end

        % Build analog1: single-sample spike at each rising edge.
        % This replaces the original O(n) element-by-element narrowing loop.
        analog1(rising) = 1;
        break;   % only one pacing channel expected
    end
end

% ── Write results into struct ─────────────────────────────────────────────────
% Dynamic field names replace the if/elseif camera_counter == 1/2/3/4 chain,
% allowing any number of cameras without code changes.
cam_f     = sprintf('CAM%d',           camera_counter);
img_f     = sprintf('CAM%d_image',     camera_counter);
raw_img_f = sprintf('CAM%d_raw_image', camera_counter);

cmos_all_data.(cam_f)     = cmosData;
cmos_all_data.(img_f)     = bgImageNom;
cmos_all_data.(raw_img_f) = bgImage;

if camera_counter == 1
    cmos_all_data.acqFreq = acqFreq;
    if ~isempty(analog1) && max(analog1) > 0
        cmos_all_data.analog1    = analog1;
        cmos_all_data.analogTime = analogTime;
    end
end
end
