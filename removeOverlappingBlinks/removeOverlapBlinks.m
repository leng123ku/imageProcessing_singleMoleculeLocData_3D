function [pixel, blinks, flag] = removeOverlapBlinks(directoryName, nimage, Tcenter, Tblinks)
%% removeOverlapBlinks removes blinks which are suspicous to be an overlapping of two or more blinks
% INPUT
% Tcenter (n_blinks,4): #file, #frame in tiff file, x_blink(px),
% y_blink(px) => later on in the program two more column will be added to
% Tcenter, actual frame number and remove/keep index.
%
%
% OUTPUT
%
%
%
%
%
%%
%fprintf('Removing potential overlapping blinks ... ');
t.removeOverlapBlinks_s = tic;
% Tcenter = Tcenter;
% Tblinks = Tblinks;

[n_blinks, n_col] = size(Tcenter);
[cut_size,~,~] = size(Tblinks);

if n_blinks ~= 0
    flag = 1;
    % formatting Tcenter2 matrix to include an absolute frame number in its fifth column.
    % currently Tcenter2(# of identified blinks, 4) =>  #tiff file, #frame, x, y
    frameNum = zeros(n_blinks,1);
    for i = 1:n_blinks
        frameNum(i,1)=sum(nimage(1:Tcenter(i,1)-1)) + Tcenter(i,2);
    end
    Tcenter(:,5) = frameNum;
    % sorting Tcenter2 according to frameNum, an ascending sorting!
    [~, order] = sort(Tcenter(:,5));
    Tcenter = Tcenter(order,:);
    Tblinks = Tblinks(:,:,order);
    % true: blink must be removed
    Tcenter(:,6) = false ;
  
    % frameInfo: frame number, number of blinks, starting row, ending row
    minFrameNum = min(Tcenter(:,5));
    maxFrameNum = max(Tcenter(:,5));
    frameInfo_noOfRow = maxFrameNum-minFrameNum+1;
    
    frameInfo = zeros(frameInfo_noOfRow,4);
    frameInfo(:,1) = (minFrameNum:maxFrameNum)';
    
    [count, ~] = hist(Tcenter(:,5), frameInfo(:,1));
    frameInfo(:,2) = count';
    frameInfo(1,3) = 1;
    frameInfo(1,4) = frameInfo(1,3)+ frameInfo(1,2)-1;
    
    for j = 2: frameInfo_noOfRow
        if frameInfo(j,2) ~= 0
            frameInfo(j,3) = sum(frameInfo(1:j-1,2)) + 1;
            frameInfo(j,4) = frameInfo(j,3) + frameInfo(j,2) - 1;
        end
    end
    % remove rows with no blinks
    %     [x,~] = find(frameInfo(:,3) ~= 0);
    %     frameInfo = frameInfo(x,:);
    
    % comparison matrix, 1st dim: #frame, 2nd dim: 1 for vicinity pixels, 2
    % for central pixel, 3rd dim:
    compPixelInfo = zeros(512,512,3);
    radius = 3;
    
    
    %     Tcenter2(9:end,:) = [];
    %     Tblinks2(:,:,9:end) = [];
 %{
    for k = 1:size(frameInfo,1)
        if frameInfo(k,2)
            for l = frameInfo(k,3):frameInfo(k,4)
                x_blink = Tcenter(l,3);
                y_blink = Tcenter(l,4);
                x_blink_neighb = (x_blink - radius) : (x_blink + radius);
                y_blink_neighb = (y_blink - radius) : (y_blink + radius);
                if k == 1
                    compPixelInfo(x_blink_neighb,y_blink_neighb,1) = frameInfo(k,1);
                    compPixelInfo(x_blink_neighb,y_blink_neighb,2) = 1;
                    %compPixelInfo(x_blink, y_blink,2) = 2;
                    compPixelInfo(x_blink-1:x_blink+1, y_blink-1:y_blink+1,2) = 2;
                    compPixelInfo(x_blink_neighb,y_blink_neighb,3) = l; % this is el, not one
                else
                    % check if there is any sequential blink, pixel ID
                    % difference = 1
                    indexCheck = frameInfo(k,1) - compPixelInfo(x_blink_neighb,y_blink_neighb,1);
                    [x,y] = find(indexCheck == 1);
                    same_blink_bool = ((compPixelInfo(x_blink,y_blink,2) == 2) && (indexCheck(radius+1, radius+1) == 1));
                    % if the blink in previous frame should be removed, so
                    % this one should be removed too.
                    l_preFrame = compPixelInfo(x_blink,y_blink,3);
                    same_blink_bool = (same_blink_bool) && ~(Tcenter(l_preFrame,6));
                    %same_blink_bool = false;
                    if ~isempty(x) && ~same_blink_bool
                        Tcenter(l,6) = true;
                        % now remove the overlapping blink in previous
                        % frame
                        for m = 1:size(x,1)
                            x_idx = x_blink + (x(m)-(radius+1));
                            y_idx = y_blink + (y(m)-(radius+1));
                            overlap_blink_pre_frame_inx = compPixelInfo(x_idx,y_idx,3);
                            % do not exclude the blinks with only one pixel
                            % difference in peak's coordinate
                            if ~((x_idx == radius+1) && (y_idx == radius+1))
                                Tcenter(overlap_blink_pre_frame_inx,6) = true;
                            end
                        end
                    end
                    compPixelInfo(x_blink_neighb,y_blink_neighb,1) = frameInfo(k,1);
                    compPixelInfo(x_blink_neighb,y_blink_neighb,2) = 1;
                    %compPixelInfo(x_blink, y_blink,2) = 2;
                    compPixelInfo(x_blink-1:x_blink+1, y_blink-1:y_blink+1,2) = 2;
                    compPixelInfo(x_blink_neighb,y_blink_neighb,3) = l;
                end
            end
        end
    end
%}
%     aa = compPixelInfo(:,:,1);
%     bb = compPixelInfo(:,:,2);
%     cc = compPixelInfo(:,:,3);
else
    error('No blink has been identified! .. Double check raw data or identification routine');
end


%% constructing the final comprehensive pixel.mat file
[x,~] = find(Tcenter(:,6) == 0);
[n_blinks_nonOverlap, ~] = size(x);
blinks_percent_removed = round( (1 -  n_blinks_nonOverlap / n_blinks) * 100);
pixel = Tcenter(x,:);
blinks = Tblinks(:,:,x);
filePath = [directoryName,'\pixel.mat'];
save(filePath,'pixel','blinks', '-v7.3')
cprintf('*red', '%g blinks condensed to %g blinks  => Overlapping blinks percentage = %g%%. \n', n_blinks, n_blinks_nonOverlap, blinks_percent_removed)
t.removeOverlapBlinks_elapsed = toc(t.removeOverlapBlinks_s);
fprintf('Elasped time for removing potentially overlapping blinks: %.0f s. \n\n', ...
    t.removeOverlapBlinks_elapsed);