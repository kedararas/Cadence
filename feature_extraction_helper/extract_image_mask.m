function [image_mask] = extract_image_mask(bgimage, iterations )

%image mask
num_rows = size(bgimage,1);
num_cols = size(bgimage,2);
mask = ones(num_rows,num_cols);
% image_bg = rgb2gray(real2rgb(bgimage, 'gray'));
image_bg = bgimage;
image_fg = activecontour(image_bg,mask,iterations);
connected_regions = bwconncomp(image_fg);
[biggest,~] = max(cellfun(@numel, connected_regions.PixelIdxList));
image_roi = bwareaopen(image_fg, round(0.25*biggest)); 
%image_roi = connected_regions.PixelIdxList{biggest};
image_mask = nan(num_rows,num_cols);
image_mask(image_roi == 1) = 1;
%image_mask = repmat(image_mask, [1 1 size(cmos_data,3)]);


end

