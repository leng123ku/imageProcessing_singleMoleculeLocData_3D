function n_blinks_tot = searchForBlinks_par_dynamicBg(n_file, directoryName, cut_xy, channel)
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
CCDgain = 25;
noiseLevel = 3;
% neigh_check: segmentation radius for blink exploration
neighCheck = 2;
neighCheckArea = 2*neighCheck+1;
% cutRadius and cutRegion for final stage, saving into center_(n_file) and
% blinks_(n_file)
cutRadius = 5;
cutRegion = 2*cutRadius + 1;
%% initiallization
maxBlinkEstimate = (21 - n_file)*1e5; % assuming maxmimum of 2 million blinks for the first file and gradually decreasing afterwards

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

n_blinks = zeros(num_image,1);
% blinkCenter_temp = zeros(200,4,num_image,'single');
% blinkCut_temp = zeros(cutRegion,cutRegion,200,num_image,'single');
blinkCenter_temp = cell(num_image,1);%zeros(maxBlinkEstimate,4,'single');
blinkCut_temp = cell(num_image,1);%zeros(cutRegion,cutRegion,maxBlinkEstimate,'single');
%% calculate the dynamic bg
cprintf('*blue','calculating the dynamic bg ... \n');
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

%sel_frames = 1:num_image/2;
x1 = cut_xy(1);
x2 = cut_xy(2);
y1 = cut_xy(3);
y2 = cut_xy(4);
dev = num_image;
n_dev = ceil( (num_image) / dev);
base = cell(num_image,1);
threshold = cell(num_image,1);

for i = 1:n_dev
    from = (i-1)*dev + 1;
    if i ~= n_dev;
        to = i*dev;
    else
        to = num_image;
    end
    sel_frames = from:to;
    n_sel_frames = size(sel_frames, 2);
    base_temp = double(zeros(width, height));
    
    %bar_base = waitbar(n_sel_frames, 'Calculating the base ...');
    parfor_progress(num_image);
    
    parfor j = from:to
        if strcmp(channel,'green')
            I =  flipud(double(imread(imageStackPath, j, 'Info', stackInfo)))/CCDgain;
        else
            I =  double(imread(imageStackPath, j, 'Info', stackInfo))/CCDgain;
        end
        I0 = I(x1:x2,y1:y2);
        base_temp = base_temp + I0;
        %waitbar(i/n_sel_frames, bar_base);
        parfor_progress;
    end
    %close(bar_base);
    parfor_progress(0);
    
    base_temp = base_temp/n_sel_frames;
    threshold_mean  = mean(base_temp(:));
    threshold_sd = std(base_temp(:));
    for k = from:to
        base{k} = base_temp;
        if strcmp(channel,'green')
            threshold{k} = 1*threshold_mean + floor(2*threshold_sd);%threshold_mean+ 2*threshold_sd; %2
            bg_sub_factor = 0;
        else
            threshold{k} = 1*threshold_mean + floor(2*threshold_sd);%threshold_mean+ 2*threshold_sd; %1.5
            bg_sub_factor = 0;
        end
    end
end


%% exploration for potential blinks
cprintf('*blue','analyzing individual frames ... \n');
parfor_progress(num_image);

parfor m = 1:2000%num_image
    if strcmp(channel,'green')
        I =  flipud(single(imread(imageStackPath, m, 'Info', stackInfo)));
    else
        I =  single(imread(imageStackPath, m, 'Info', stackInfo));
    end

    A0 = I(x1:x2,y1:y2)/CCDgain;
    %A0 = I(cut_xy(1):cut_xy(2), cut_xy(3):cut_xy(4))/CCDgain;
    %base = 0.9*base + 0.099*min(base, A0) + 0.001*A0;
    
    
    A = A0 - bg_sub_factor*base{m};
    
    %A = A0 - threshold{m};
    %A = medfilt2(A, [2 2]);
    %A = A0 - base;
    A( A<0 ) = 0;
    
    %A = wiener2( A, [5 5] ); % previously using [5 5]
    %A = imfilter(A,h0,'replicate');
    %[A, Abg] = imag_background_sub(A1,h2);
    flag = false(width, height);
    
    flag(A>threshold{m}) = 1;
    
    flag(1:2*cutRadius,:)= false;
    flag( width-2*cutRadius : width, :)= false;
    flag(:,1:2*cutRadius)= false;
    flag(:, height-2*cutRadius : height)= false;
    
    flag_sum = sum(flag(:));
    
    if flag_sum >= 1
        %% calculate maxima for every (neighCheckArea * neighCheckArea) image to identify the center
        [x_flag,y_flag]= find(flag > 0);
        [num_bead_beforeMaxCheck,~]= size(x_flag); %flag_sum = num_bead here
        for i = 1:num_bead_beforeMaxCheck
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
        [num_bead_afterMaxCheck,~]= size(x_flag);

        if (num_bead_afterMaxCheck >= 2)
            for i=1:num_bead_afterMaxCheck-1
                if flag(x_flag(i),y_flag(i))
                    for j=i+1:num_bead_afterMaxCheck
                        if flag(x_flag(j),y_flag(j))
                            dist=sqrt((x_flag(i)-x_flag(j))^2+(y_flag(i)-y_flag(j))^2);
                            if (dist <= neighCheck+3)
                                if A(x_flag(i),y_flag(i)) > A(x_flag(j),y_flag(j))
                                    flag(x_flag(j),y_flag(j))=  false;
                                else
                                    flag(x_flag(i),y_flag(i))= false;
                                end
                            end
                        end
                    end
                end
            end
        end
        
        [x_flag,y_flag] = find(flag > 0);
        [num_bead_afterNeighbRemoval,~] = size(x_flag);  % final means after neighbors removal
        if (num_bead_afterNeighbRemoval >= 1)
            for i= 1:num_bead_afterNeighbRemoval
                
                I = A0(x_flag(i)-cutRadius:x_flag(i)+cutRadius, y_flag(i)-cutRadius:y_flag(i)+cutRadius);
                %I = imgaussfilt(I, cutRadius);
                [~,~,~,wx0,wy0] = gauss2dellipse(I,xMesh,yMesh, noiseLevel);
                %wx0 = 2; wy0 = 2;
                
                if (wx0 > 0.4) && (wy0 > 0.4) % 1 1
                    if (wx0 < 12) && (wy0 < 12)  %8 8
                        
                        n_blinks(m) = n_blinks(m) + 1;
                        % blinkCenter => n_file, frame#, x(px), y(px)
                        blink_id = [n_file m x_flag(i) y_flag(i)];
                        % blinkCut: cuts around identified beads
                        if i == 1
                            blinkCenter_temp{m} = {blink_id};
                            blinkCut_temp{m} = {I};
                        else
                            blinkCenter_temp{m} = [blinkCenter_temp{m}; {blink_id}];
                            blinkCut_temp{m} = [blinkCut_temp{m}; {I}];
                        end
                    end
                else
                    flag(x_flag(i), y_flag(i)) = 0;
                end
            end
        end   
    end % flag_sum >= 1 
    parfor_progress;
end
parfor_progress(0);

% figure; imshow(flag, [0 1]);
% figure; imshow(A0, [min(A0(:)) max(A0(:))]);
% figure; imshow(A, [min(A(:)) max(A(:))]);
%%
% n_blinks = sum(n_blinks);
% blinkCenter = zeros(n_blinks,4,'single');
% blinkCut = zeros(cutRegion,cutRegion,n_blinks,'single');
% 
% n = 0;
% for i = 1:num_image
%     for j = 1:200
%         if sum(any(blinkCenter_temp(j,:,i)))
%             n = n+1;
%             blinkCenter(n,:) = blinkCenter_temp(j,:,i);
%             blinkCut(:,:,n) = blinkCut_temp(:,:,j,i);
%         end 
%     end
% end

%% accumulate results and save
n_blinks_tot = sum(n_blinks);
blinkCenter = zeros(n_blinks_tot,4,'single');
blinkCut = zeros(cutRegion,cutRegion,n_blinks_tot,'single');
idx = 0;

for i = 1:num_image
    for j = 1:n_blinks(i)
        idx = idx+1;
        blinkCenter(idx,:) = blinkCenter_temp{i}{j};
        blinkCut(:,:,idx) = blinkCut_temp{i}{j};
    end
end

center = blinkCenter(1:n_blinks_tot,:);
blinks = blinkCut(:,:,1:n_blinks_tot);
clear blinkCenter blinkCut
filePath_pixel_n = [directoryName,'\center',int2str(n_file),'.mat'];
save(filePath_pixel_n,'center','blinks');
t1 = toc(t0); % elapsed time in min.
fps = num_image / t1; % frame/sec
cprintf('*String', ' %12s  =>  %7g   blinks found in %g frames in %2.1f min   =>   %2.0f FPS. \n', fprintf_idx, n_blinks_tot, num_image, t1/60, fps);
