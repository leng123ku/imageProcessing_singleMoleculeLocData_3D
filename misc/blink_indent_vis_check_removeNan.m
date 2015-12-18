function blink_indent_vis_check_removeNan(directoryName, sel_file, sel_frame, channel)

imageStackPath = fullfile(directoryName,'STORM.tif');
stackInfo = imfinfo(imageStackPath);
%num_image = numel(stackInfo);
if strcmp(channel,'green')
    I =  flipud(double(imread(imageStackPath, sel_frame, 'Info', stackInfo)));
else
    I =  double(imread(imageStackPath, sel_frame, 'Info', stackInfo));
end
filePath = fullfile(directoryName, 'position.mat');
load(filePath); % position.mat
rows_to_catch = ones(size(position,1),1);
rows_to_catch((position(:,14) ~= sel_file) | (position(:,15) ~= sel_frame)) = NaN;
position(isnan(rows_to_catch),:) = [];
iden_pt = position(:,[3 5]); % x(nm) y(nm)

fig1 = figure; set(fig1,'name', 'fluor_check_removeNans'); imshow(I,[1e2 0.5*max(I(:))]); colormap hot;
viscircles([iden_pt(:,2) iden_pt(:,1)], 10*ones(size(iden_pt,1),1), 'edgecolor','b');
pause(1);
%waitforbuttonpress;