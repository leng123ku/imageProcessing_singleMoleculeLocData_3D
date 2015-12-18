function blink_indent_vis_check(directoryName, sel_file, sel_frame, Tcenter, channel)

imageStackPath = fullfile(directoryName,'STORM.tif');
stackInfo = imfinfo(imageStackPath);
%num_image = numel(stackInfo);
if strcmp(channel,'green')
    I =  flipud(double(imread(imageStackPath, sel_frame, 'Info', stackInfo)));
else
    I =  double(imread(imageStackPath, sel_frame, 'Info', stackInfo));
end
%I =  double(imread(imageStackPath, sel_frame, 'Info', stackInfo));
rows_to_catch = ones(size(Tcenter,1),1);
rows_to_catch((Tcenter(:,1) ~= sel_file) | (Tcenter(:,2) ~= sel_frame)) = NaN;
Tcenter(isnan(rows_to_catch),:) = [];
iden_pt = Tcenter(:,[3 4]); % x(nm) y(nm)

fig1 = figure; set(fig1,'name', 'fluor_check'); imshow(I,[1e2 1*max(I(:))]); colormap hot; %2e3
viscircles([iden_pt(:,2) iden_pt(:,1)], 10*ones(size(iden_pt,1),1), 'edgecolor','b');
pause(1);
%waitforbuttonpress;