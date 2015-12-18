function [n_file, n_image, imageSize_x, imageSize_y]  = searchForSMImageStacks(directoryName)
%% searchForSMImageStacks finds the available .tif images stacks in directoryName path
% SM image stacks must be named as STORM.tif, STORM_X2.tif, ...
%{ 
	- INPUTS:
        - directoryName: the path to SM image stacks
    - OUTPUTS:
        - n_file: number of .tif files found in the path
        - n_image: a row vector whose elements correspond to number of
        frames in each image stack.
        - imageSize_x and imageSize_y: frame size in pixel
%}
filePath=[directoryName,'\STORM','.tif'];
flag1=exist(filePath, 'file');
n_image = zeros(1,20); 
if flag1 == 0
    error('First SM image stack is supposed to be named as STORM.tif => %s does not exist. \n',filePath);
else
    n_file = 1;
    I = (filePath);
    info = imfinfo(I);
    n_image(1) = numel(info);
    data = imread(I, 1, 'Info', info);
    [imageSize_x, imageSize_y] = size(data);
end

for nf = 2:30
    filePath = [directoryName,'\STORM_X',int2str(nf),'.tif'];
    flag1 = exist(filePath,'file');
    if flag1 ~= 0
        n_file = nf;
        I = (filePath);
        info = imfinfo(I);
        n_image(nf) = numel(info);
    end
end
% removing extra columns in nimage
n_image(:, n_file+1:20) = [];