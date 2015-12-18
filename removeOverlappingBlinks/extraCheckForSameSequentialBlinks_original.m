function extraCheckForSameSequentialBlinks(directoryName, nimage, CCDpixelx,CCDpixely)

t.removeOverlapBlinks_s = tic;
filePath= fullfile(directoryName,'\position_fit.mat');
load(filePath);
n_blinks = size(position_fit,1);
% xyFrame: x y #file #frameInFile
xyFrame = zeros(n_blinks,6);
xyFrame(:,1:4) = [position_fit(:,3)*CCDpixelx position_fit(:,5)*CCDpixely position_fit(:,14) position_fit(:,15)];
for i = 1:n_blinks
    xyFrame(i,5)=sum(nimage(1:position_fit(i,14)-1)) + position_fit(i,15);
end
% true: blink must be removed
xyFrame(:,6) = false;

% frameInfo: frame number, number of blinks, starting row, ending row
minFrameNum = min(xyFrame(:,5));
maxFrameNum = max(xyFrame(:,5));
frameInfo_noOfRow = maxFrameNum-minFrameNum+1;

frameInfo = zeros(frameInfo_noOfRow,4);
frameInfo(:,1) = (minFrameNum:maxFrameNum)';

[count, ~] = hist(xyFrame(:,5), frameInfo(:,1));
frameInfo(:,2) = count';
frameInfo(1,3) = 1;
frameInfo(1,4) = frameInfo(1,3)+frameInfo(1,2)-1;

for j = 2: frameInfo_noOfRow
    if frameInfo(j,2) ~= 0
        frameInfo(j,3) = sum(frameInfo(1:j-1,2)) + 1;
        frameInfo(j,4) = frameInfo(j,3) + frameInfo(j,2) - 1;
    end
end

%%{
for i = 2:n_blinks-1
    if (xyFrame(i,5) > minFrameNum+1) && (xyFrame(i,5) < maxFrameNum-1)
        % check if there is any event in previous or next frame
        [frameInx, ~] = find(frameInfo(:,1) == xyFrame(i,5));
        framePre_bool = frameInfo(frameInx-1,2) ~= 0;
        frameNex_bool = frameInfo(frameInx+1,2) ~= 0;
%         framePrePre_bool = frameInfo(frameInx-2,2) ~= 0;
%         frameNexNex_bool = frameInfo(frameInx+2,2) ~= 0;
        
        if (framePre_bool || frameNex_bool)
            if framePre_bool
                dxy = pdist2(xyFrame(i,[1 2]), xyFrame(frameInfo(frameInx-1,3):frameInfo(frameInx-1,4),[1 2]));
                [distThresholdInx, ~] = find(dxy(:) < 400 & dxy(:) > 40);
                if isempty(distThresholdInx)
                    framePre_bool = false;
                end
            end
            
            
            if frameNex_bool
                dxy = pdist2(xyFrame(i,[1 2]), xyFrame(frameInfo(frameInx+1,3):frameInfo(frameInx+1,4),[1 2]));
                [distThresholdInx, ~] = find(dxy < 400 & dxy > 40);
                if isempty(distThresholdInx)
                    frameNex_bool = false;
                end
            end
%{          
            if framePrePre_bool
                dxy = pdist2(xyFrame(i,[1 2]), xyFrame(frameInfo(frameInx-2,3):frameInfo(frameInx-2,4),[1 2]));
                [distThresholdInx, ~] = find(dxy(:) < 400 & dxy(:) > 50);
                if isempty(distThresholdInx)
                    framePrePre_bool = false;
                end
            end
            
            if frameNexNex_bool
                dxy = pdist2(xyFrame(i,[1 2]), xyFrame(frameInfo(frameInx+2,3):frameInfo(frameInx+2,4),[1 2]));
                [distThresholdInx, ~] = find(dxy < 400 & dxy > 50);
                if isempty(distThresholdInx)
                    frameNexNex_bool = false;
                end
            end
%}          
            
            
            %if (framePre_bool || frameNex_bool) || (framePrePre_bool || frameNexNex_bool)
             if (framePre_bool || frameNex_bool) 
                xyFrame(i,6) = true;
                %xyFrame(frameInfo(frameInx+1,3) + distThresholdInx,6) = true;
                %xyFrame(frameInfo(frameInx-1,3) + distThresholdInx,6) = true;
            end
        end
    end
end
%%}

[rowsToKeep, ~] = find(xyFrame(:,6) == false);
[n_blinks_nonOverlap, ~] = size(rowsToKeep);
position_fit_temp = position_fit(rowsToKeep,:);
% converting the accumulative total number of blinks to an averaged
% amplitude => A/(Wx*Wy)
position_fit_temp (:,2) = position_fit_temp(:,2)./( position_fit_temp(:,4).* position_fit_temp(:,6));
position_fit_falseSequentialBlinksRemoved = position_fit_temp;
filePath_position_fit = [directoryName,'\position_fit_falseSequentialBlinksRemoved.mat'];
save(filePath_position_fit , 'position_fit_falseSequentialBlinksRemoved');
cprintf('*red', '%g blinks condensed to %g blinks. \n',n_blinks, n_blinks_nonOverlap)
t.removeOverlapBlinks_elapsed = toc(t.removeOverlapBlinks_s);
fprintf('Elasped time for removing potentially overlapping blinks: %.0fs. \n\n', ...
    t.removeOverlapBlinks_elapsed);
