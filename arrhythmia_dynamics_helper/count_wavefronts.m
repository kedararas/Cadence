function [WFfront,WFcount] = count_wavefronts(phase_data, wf_size_threshold, pixel_size_mm, min_wf_length_mm)
% count_wavefronts  Detect and count wavefronts in phase map data.
%
%   [WFfront, WFcount] = count_wavefronts(phase_data, wf_size_threshold)
%       Uses wf_size_threshold directly (pixels) as the minimum contour
%       length. Falls back to 10% of frame diagonal if empty/zero.
%
%   [WFfront, WFcount] = count_wavefronts(phase_data, [], pixel_size_mm)
%       Derives threshold from pixel size (mm/pixel) and a default minimum
%       wavefront length of 3 mm.
%
%   [WFfront, WFcount] = count_wavefronts(phase_data, [], pixel_size_mm, min_wf_length_mm)
%       Full physical specification: threshold = min_wf_length_mm / pixel_size_mm.

%% Resolve wavefront size threshold
frame_diagonal = sqrt(size(phase_data,1)^2 + size(phase_data,2)^2);

if nargin >= 3 && ~isempty(pixel_size_mm) && pixel_size_mm > 0
    % Physical specification: derive threshold from real-world units
    if nargin < 4 || isempty(min_wf_length_mm)
        min_wf_length_mm = 3; % default: reject wavefronts shorter than 3 mm
    end
    wf_size_threshold = round(min_wf_length_mm / pixel_size_mm);
    fprintf('count_wavefronts: pixel_size=%.3f mm, min_wf=%.1f mm -> threshold=%d px\n', ...
        pixel_size_mm, min_wf_length_mm, wf_size_threshold);
elseif nargin < 2 || isempty(wf_size_threshold) || wf_size_threshold <= 0
    % Fallback: 15% of frame diagonal
    wf_size_threshold = round(0.15 * frame_diagonal);
    fprintf('count_wavefronts: no threshold supplied, using 15%% of diagonal -> %d px\n', ...
        wf_size_threshold);
end

%% Calculate Wavefront***
level=pi/2; % wavefront definition
num_frames = size(phase_data, 3);

% Pre-allocate outputs to avoid repeated reallocation inside the loop
WFfront = cell(2, num_frames);
WFcount = zeros(2, num_frames);

for w=1:num_frames % loop through time
    % Reset per-frame accumulators at top of loop
    Wave = {};
    l = 1;

    %Step 1: Use contour algorithm to calculate wavefronts
    testWave = contourcs(phase_data(:,:,w), [level level]);

    %Step 2: Remove any wavefronts <10 pixels in length - noise filter
    threshWave = wf_size_threshold;  % change this value to modify wavefront visualization
    num_contours = numel(testWave);

    for c = 1:num_contours
        seg = testWave(c);          % extract struct once to avoid repeated indexing
        if seg.Length > threshWave
            Wave{l} = [seg.X; seg.Y];
            l = l + 1;
        end
    end

    %Step 3: Remove and separate wavefronts based on large spatial gradient
    [px, py] = gradient(phase_data(:,:,w));
    spaceGrad = max(abs(px), abs(py));
    [testerx, testery] = find(spaceGrad > 1);
    test = [testery, testerx];
    segWave = {};
    index = 1;
    for r = 1:numel(Wave) % loop through all waves in contour structure
        X = Wave{r};
        tempX_floor = floor(X'); % test both floor and ceiling to deal with rounding errors
        tempX_ceil  = ceil(X');
        % Combine both into one ismember call to halve the search cost
        testcomp_both = ismember([tempX_floor; tempX_ceil], test, 'rows');
        testcomp_floor = testcomp_both(1:end/2);
        testcomp_ceil  = testcomp_both(end/2+1:end);
        testcomp = max(testcomp_ceil, testcomp_floor)';
        Y = X(:, ~testcomp);

        %Separate non connected vectors
        if size(Y, 2) < 2
            continue
        end
        Yshift = Y(:, 2:end);
        dist = sqrt((Yshift(1,:) - Y(1,1:end-1)).^2 + (Yshift(2,:) - Y(2,1:end-1)).^2);
        segmentpts = find(dist > 5);
        if isempty(segmentpts)
            segWave{index} = Y;
            index = index + 1;
        else
            segWave{index} = Y(:, 1:segmentpts(1)); index = index + 1;
            if length(segmentpts) > 1
                for ii = 1:length(segmentpts)-1
                    segWave{index} = Y(:, segmentpts(ii)+1:segmentpts(ii+1));
                    index = index + 1;
                end
            end
            segWave{index} = Y(:, segmentpts(end)+1:end);
            index = index + 1;
        end

    end

    num_wf = numel(segWave);
    wf_count = 0;
    newsegWave = {};
    modified_wf_size = zeros(1, num_wf);
    for i = 1:num_wf
        if ~isempty(segWave{i}) && size(segWave{i}, 2) >= 2*wf_size_threshold
            modified_wf_size(i) = size(segWave{i}, 2);
            wf_count = wf_count + 1;
            wf = segWave{i};
            if wf(1,1) == wf(1,end) && wf(2,1) == wf(2,end)
                newsegWave{1,wf_count} = wf(:, 1:end-1);
            else
                newsegWave{1,wf_count} = segWave{i};
            end
        end
    end

    % Filter empties without cellfun by only keeping non-empty entries
    nonempty_mask = ~cellfun('isempty', segWave);
    WFfront{1,w} = newsegWave;
    WFfront{2,w} = segWave(nonempty_mask);
    WFcount(1,w) = sum(modified_wf_size >= 4*wf_size_threshold);
    WFcount(2,w) = wf_count;

end

% Plot results
% if video_mode
%     fig = figure;
%     video_file = cmos_all_data.file_name;
%     video_file = strrep(video_file, '.mat', '-phase.mp4');
%     video_file = strrep(video_file, '/arrhythmia/','/arrhythmia/movies/');
%     writerObj = VideoWriter(video_file, 'MPEG-4');
%     writerObj.FrameRate = 20; % change this value to modify speed of movie
%     open(writerObj);
%     movegui(fig,'center');
%     set(gcf, 'color', [1 1 1]);
% 
%     image_file = strrep(video_file, '.mp4', '.png');
%     image_file = strrep(image_file, '/movies/','/images/phase/');
% 
%     frame_size = size(phase_data,3);
%     if frame_size > 4000
%         frame_size = 4000;
%     end
%     for i = 1:1:frame_size
%         G = cmos_all_data.bgimage;
%         image(G);
%         hold on;
%         Mframe = phase_data(:,:,i);
%         contourf(Mframe, 32, 'LineStyle', 'none');
%         colormap jet;
%         %contourcbar;
%         caxis([-pi pi])
%         %axis image
%         axis off
%         hold on
% 
%         if exist('segWave','var')
%             segWave=WFfront{1,i};
%             for pp=1:size(segWave,2)
%                 WF=segWave{1,pp};
%                 plot(WF(1,:), WF(2,:), 'w', 'LineWidth', 2);
%             end
%         end
%         pause(.001)
%         hold off
%         frame = getframe(gcf);
%         writeVideo(writerObj,frame);
% 
% %         if i < 10
% %             prefix = strcat('-000',num2str(i), '.png');
% %         elseif i < 100
% %             prefix = strcat('-00',num2str(i), '.png');
% %         elseif i < 1000
% %             prefix = strcat('-0',num2str(i), '.png');
% %         else
% %             prefix = strcat('-',num2str(i),'.png');
% %         end
% %         test_file = strrep(image_file, '.png', prefix);
% %         saveas(fig, test_file);
% 
%     end
%     close(fig);
%     close(writerObj);
%     close all;
% end
% 
% 
% %Step 3: Generate Phase Activation Map (activation is pi/2)
% %     num_frames = size(phase_data,3);
% %     phase_act_data = nan(100,100,num_frames);
% %     for i = 1:size(phase_data,3)
% %         phase_frame = round(squeeze(phase_data(:,:,i)).* 10)/10;
% %         phase_mask = phase_frame;
% %         phase_mask(~isnan(phase_mask)) = 1;
% %         phase_frame(phase_frame ~= 1.6) = 0;
% %         phase_frame(phase_frame == 1.6) = 1;
% %         phase_frame = phase_frame .* phase_mask;
% %         phase_act_data(:,:,i) = phase_frame;
% %     end
% 
% message = strcat('Extracted Wavefront Count');
% disp(message);
% 
% end

