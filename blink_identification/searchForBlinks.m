function n_blinks = searchForBlinks(n_file, directoryName, cut_xy)
%% searchForBlinks analyzes .tiff file named "STORM_X(n_file)" to find potential blinks
% results are saved into two .mat files
% center_(n_file) => (n_file, frame_number, x(px), y(px))
% blinks_(n_file) => (cutRegion, cutRegion, n_blinks)     segments of starting image stack containing potential blinks
% 3rd dimension of blinks corresponds to row# of center.
%{
INPUTS:
n_file: file number to analyze
directoryName: files' path
cut_xy =>  cut_xy(1:4) : cutting area to analyze; top, bottom, left and
right pixels respectively.

%}
%% user's settings
t0 = tic;
threshold = 10;%12
CCDgain = 250;
noiseLevel = 2;
% neigh_check: segmentation radius for blink exploration
neighCheck = 2;
neighCheckArea = 2*neighCheck+1;
% cutRadius and cutRegion for final stage, saving into center_(n_file) and
% blinks_(n_file)
cutRadius = 5;
cutRegion = 2*cutRadius + 1;
%% initiallization
maxBlinkEstimate = (21 - n_file)*1e5; % assuming maxmimum of 2 million blinks for the first file
blinkCenter = zeros(maxBlinkEstimate,4,'single');
blinkCut = zeros(cutRegion,cutRegion,maxBlinkEstimate,'single');

[yMesh, xMesh] = meshgrid(1:cutRegion);
% mesh(:,1) = xMesh(:);
% mesh(:,2) = yMesh(:);

% finiding the image size
width = cut_xy(2) - cut_xy(1) + 1;
height = cut_xy(4) - cut_xy(3) + 1;

base = ones(width, height, 'single') * 1e7;
h0 = ones(2,2)/4; % smoothing mask
n_blinks = 0;

if (n_file == 1)
    imageStackPath = fullfile(directoryName,'STORM.tif');
    fprintf_idx = 'STORM';
else
    imageStackPath =[directoryName '\STORM_X' int2str(n_file) '.tif'];
    fprintf_idx = ['STORM_X' int2str(n_file)];
end

stackInfo = imfinfo(imageStackPath);
num_image = numel(stackInfo);

%% calculate a baseline image using the first 300 image of each stack
% for i = 1:2:301
%     I =  single(imread(imageStackPath, i, 'Info', stackInfo));
%     I_cut = I(cut_xy(1):cut_xy(2), cut_xy(3):cut_xy(4))/CCDgain;
%     base = min(base, I_cut);
% end
% % another iteration
% for i = 1:3
%     for j = 1:301
%         I =  single(imread(imageStackPath, j, 'Info', stackInfo));
%         I_cut = I(cut_xy(1):cut_xy(2), cut_xy(3):cut_xy(4))/CCDgain;
%         base = 0.9*base + 0.099*min(base, I_cut) + 0.001*I_cut;
%     end
% end

base = double(zeros(width, height));
sel_frames = 1:4:num_image;
n_sel_frames = size(sel_frames, 2);
%     for i = frames_to_probe_base
%         I =  double(imread(info.imageStackPath, i, 'Info', info.stackInfo))/info.CCDgain;
%         base = base + I;
%     end

for i = 1:n_sel_frames
    I =  double(imread(imageStackPath, sel_frames(i), 'Info', stackInfo))/CCDgain;
    I0 = I(cut_xy(1):cut_xy(2), cut_xy(3):cut_xy(4));
    base = base + I0;
end

%base = base / (NumFramBead_Down + NumFramBead_Up + 1);
base = base / n_sel_frames;
threshold_sd = std(base(:));
threshold_mean = mean(base(:));
threshold = 20;%threshold_mean;% + 2*threshold_sd;

%% explorarion for potential blinks
for m = 1:num_image
    I =  single(imread(imageStackPath, m, 'Info', stackInfo));
    A0 = I(cut_xy(1):cut_xy(2), cut_xy(3):cut_xy(4))/CCDgain;
    %base = 0.9*base + 0.099*min(base, A0) + 0.001*A0;
    A = A0 - base;
    %A = A0;
    %A = medfilt2(A, [2 2]);
    %A = A0 - base;
    A( A<0 ) = 0;
    %A = wiener2( A, [5 5] ); % previously using [5 5]
    %A = imfilter(A,h0,'replicate');
    %[A, Abg] = imag_background_sub(A1,h2);
    flag = false(width, height);
    flag( A>threshold ) = 1;
    
    flag(1:2*cutRadius,:)= false;
    flag( width-2*cutRadius : width, :)= false;
    flag(:,1:2*cutRadius)= false;
    flag(:, height- 2*cutRadius : height)= false;
    
    flag_sum = sum(flag(:));
    
    if flag_sum >= 1
        %% calculate maxima for every (neighCheckArea * neighCheckArea) image to identify the center
        [x_flag,y_flag]= find(flag > 0);
        [num_bead.beforeMaxCheck,~]= size(x_flag); %flag_sum = num_bead here
        for i = 1:num_bead.beforeMaxCheck
            n1x = x_flag(i)-neighCheck;
            n2x = x_flag(i)+neighCheck;
            n1y = y_flag(i)-neighCheck;
            n2y = y_flag(i)+neighCheck;
            
            %             B = A(n1x:n2x,n1y:n2y);
            %             B (B < max(B(:))) = false;
            %             flag(n1x:n2x,n1y:n2y) = B;
            B = A(n1x:n2x,n1y:n2y);
            B_sort = sort(B(:));
            B (B < B_sort(end)) = false;
            flag(n1x:n2x,n1y:n2y) = B;
            
            %A_max = max(max(A(n1x:n2x,n1y:n2y)));
            %B_max = max(B(:)); % this is several time faster compared to the line above!
            %[ind_x, ind_y] = find(B == B_max);
            %flag(n1x:n2x,n1y:n2y) = 0;
            %flag(ind_x+n1x-1, ind_y+n1y-1) = 1;
        end
        %clear ix n1x n2x iy n1y n2y A_max i
        % Zero edges, necessary as some zeros may be pickes as max in above
        % loop!
        flag(1:2*cutRadius,:)= false;
        flag( width-2*cutRadius : width, :)= false;
        flag(:,1:2*cutRadius)= false;
        flag(:, height- 2*cutRadius : height)= false;
        %% neighbour removal: treat blinks within a pixel as the same blink
        %clear x_flag y_flag
        [x_flag,y_flag]= find(flag > 0);
        [num_bead.afterMaxCheck,~]= size(x_flag);
        if (num_bead.afterMaxCheck >= 2)
            for i=1:num_bead.afterMaxCheck-1
                for j=i+1:num_bead.afterMaxCheck
                    dist=sqrt((x_flag(i)-x_flag(j))^2+(y_flag(i)-y_flag(j))^2);
                    if (dist <= neighCheck+1)
                        if A(x_flag(i),y_flag(i)) > A(x_flag(j),y_flag(j))
                            flag(x_flag(j),y_flag(j))=  false;
                        else
                            flag(x_flag(i),y_flag(i))= false;
                        end
                    end
                end
            end
        end
        
        [x_flag,y_flag] = find(flag > 0);
        [num_bead.afterNeighbRemoval,~] = size(x_flag);  % final means after neighbors removal
        if (num_bead.afterNeighbRemoval >= 1)
            for i= 1:num_bead.afterNeighbRemoval
                
                I = A(x_flag(i)-cutRadius:x_flag(i)+cutRadius, y_flag(i)-cutRadius:y_flag(i)+cutRadius);
                [~,~,~,wx0,wy0] = gauss2dellipse(I,xMesh,yMesh, noiseLevel);
                
                if (wx0 > 0.8) && (wy0 > 0.8)
                    if (wx0 < 6) && (wy0 < 6)
                        if ~((wx0 > 5) && (wy0 > 5))
                            
                            n_blinks = n_blinks + 1;
                            % blinkCenter => n_file, frame#, x(px), y(px)
                            blinkCenter (n_blinks, 1) = n_file;
                            blinkCenter (n_blinks, 2) = m; % #frame in the current image stack
                            blinkCenter (n_blinks, 3) = x_flag(i);
                            blinkCenter (n_blinks, 4) = y_flag(i);
                            % blinkCut: cuts around identified beads
                            blinkCut(:,:,n_blinks) = I;
                            
                        end
                    end
                end
            end
        end   
    end % flag_sum >= 1 
end

%% accumulate results and save
center = blinkCenter(1:n_blinks,:);
blinks = blinkCut(:,:,1:n_blinks);
clear blinkCenter blinkCut
filePath_pixel_n = [directoryName,'\center',int2str(n_file),'.mat'];
save(filePath_pixel_n,'center','blinks');
t1 = toc(t0); % elapsed time in min.
fps = num_image / t1; % frame/sec
cprintf('*String', ' %12s  =>  %7g   blinks found in %g frames in %2.1f min   =>   %2.0f FPS. \n', fprintf_idx, n_blinks, num_image, t1/60, fps);
