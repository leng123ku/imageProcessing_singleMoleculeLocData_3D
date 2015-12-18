function [Tcenter, Tblinks] = combine_blinks(n_file, n_blinks_total, directoryName)
%% combine_blinks combines all the blinks found and stored in separate Tcenter_n and Tblinks_n .m files

filePath_pixel_n = [directoryName,'\center1.mat'];
load(filePath_pixel_n);
[n_blinks1, np] = size(center);
[n_pixelx, n_pixely, n_blinks2] = size(blinks);
if n_blinks1 ~= n_blinks2
    error('Mismatch: center.m and blinks.m => number of blinks do not match');
end
Tcenter = zeros(n_blinks_total, np, 'single');
Tblinks = zeros(n_pixelx, n_pixely, n_blinks_total, 'single');
%n1=1;
n2=0;
for nf = 1:n_file
    filePath_pixel_n = [directoryName,'\center',int2str(nf),'.mat'];
    load(filePath_pixel_n,'center','blinks')
    [n_blinks, np]=size(center);
    n1 = n2+1;
    n2 = n2+ n_blinks;
    Tcenter(n1:n2,:) = center;
    Tblinks(:,:,n1:n2) = blinks;
    clear center blinks
    delete(filePath_pixel_n)
end
% filePath_pixel_n = [directoryName,'\Tcenter.mat'];
% save(filePath_pixel_n,'Tcenter','Tblinks')