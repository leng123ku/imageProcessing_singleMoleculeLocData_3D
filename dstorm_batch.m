%% dstorm: an script to analyze dstorm data

clc; clear all; close all; close all hidden; diary ON;
warning('off','MATLAB:rankDeficientMatrix');
set(0,'DefaultFigureWindowStyle','docked');

info;
num_setOfdate = size(date, 1);
cprintf('*blue', 'Checking whether or not all required files EXIST ...\n');
check_required_files_exist;

for k = 1:num_setOfdate
    this_date = date(k,:);
    this_cellDescrip = cellDescrip(k,:);
    this_runNum_temp = runNum(k,:);
    this_runNum = this_runNum_temp (this_runNum_temp ~= 0);

    %% user's setting => creating s(setParam) object
    % setParam(date, cellDescrip, runNUm, compOrPlot)
    s = setParam(this_date, this_cellDescrip, this_runNum, 1);
    
    cprintf('*blue', 'Copying the calibration files into %s \n', s.imageStacksPathRoot);
    s.creat_copy_calibFiles(); fprintf('\n');
    
    %% startParallelPool(n_workers)
    cprintf('*blue', 'Starting a parallel pool with %g workers ...\n', s.n_workers);
    startParallelPool(s.n_workers); fprintf('\n');
    
    %% analysis
    for i = 1:numel(s.runNum)
        t_dstorm = tic;
        cprintf('*red', 'Analyzing run#%u => %s \n\n', s.runNum(i), s.imageStacksPathByRuns{i});
        
        % searchForSMImageStacks
        cprintf('*blue', 'Searching for SM image stacks ...\n');
        [n_file, n_image, imageSize_x, imageSize_y]  = searchForSMImageStacks(s.imageStacksPathByRuns{i});
        cprintf('*String','%g image stacks found. \n', n_file);
        
        % determine image size and cut off only the cell region
        cut_xy  = cut_image(s.imageStacksPathByRuns{i}, imageSize_x, imageSize_y, 0); % 0 is the cut_edge
        %cut_xy = [1,512,1,512];
        %cut_xy = [1,30,1,30];
        fprintf('Images cut at (%g, %g) and (%g, %g) pixels. \n\n', cut_xy);
        %
        if s.compOrPlot == 1
            cprintf('*blue', 'Searching for beads used for tracking + drift calculations ... \n');
            n_beads = searchForBeads(s.imageStacksPathByRuns{i});
            cprintf('*String','%g beads used for drift correction. \n', n_beads);
            
            % fit beads to calculate drift
            [Drift, nframe, bead_xyz]  = drift_bead(s.imageStacksPathByRuns{i}, ...
                n_beads, s.calFileBead, s.CCDpixely_SM, s.CCDpixelx_SM, ...
                s.CCDpixely_beadCam, s.CCDpixely_beadCam, n_image);
            cprintf('*String','Drift calculation DONE. \n\n');
            
            % searchForBlinks
            cprintf('*blue','Searching for potential blinks => Image processing ... \n');
            n_blinks_total = 0;
            for j = 1:n_file
                %n_blinks = searchForBlinks_par(j, s.imageStacksPathByRuns{i}, cut_xy);
                n_blinks = searchForBlinks_par_dynamicBg(j, s.imageStacksPathByRuns{i}, cut_xy, channel);
                n_blinks_total =  n_blinks_total +  n_blinks;
            end
            cprintf('*red', '%g potential blinks were found in %g images. \n', n_blinks_total, sum(n_image));
            
            % accumulates blinks found
            [Tcenter, Tblinks] = combine_blinks(n_file, n_blinks_total, s.imageStacksPathByRuns{i});
            cprintf('*String', 'Blinks accumulation DONE. \n\n');
            % Visual check_blink identification
            blink_indent_vis_check(s.imageStacksPathByRuns{i}, 1, 1000, Tcenter, channel);
            
            % remove overlapping blinks
            cprintf('*blue', 'Removing overlapping blinks ... \n');
            [pixel, blinks, flag_blinks]  = removeOverlapBlinks(s.imageStacksPathByRuns{i}, n_image, Tcenter, Tblinks);
            
            % fitting the blinks
            cprintf('*blue', 'Fitting the blinks ... \n');
            [~] = fit_blinks(s.imageStacksPathByRuns{i});
            blink_indent_vis_check_fitting(s.imageStacksPathByRuns{i}, 1, 1000, channel);
            % Extra check to remove same sequential blinks
            cprintf('*blue', 'Extra check to remove identical sequential blinks ... \n');
            extraCheckForSameSequentialBlinks(s.imageStacksPathByRuns{i}, n_image, s.CCDpixelx_SM, s.CCDpixely_SM);
            blink_indent_vis_check_sequential_removed(s.imageStacksPathByRuns{i}, 1, 1000, channel);
            % Remove bad fittings (NaN) and add drift corrections
            cprintf('*blue', 'Remove blinks with any NaN fitting parameter + add drift corrections ... \n');
            [~] = Remove_bad_fitting(s.imageStacksPathByRuns{i}, n_file, n_image);
            blink_indent_vis_check_removeNan(s.imageStacksPathByRuns{i}, 1, 1000, channel);
        end
        
        % Z locallization
        cprintf('*blue', 'Z-localization, filtering, plotting and writing the ViSP file \n\n');
        cprintf('*String', 'Using %s callibration curve. \n' , s.calFileSMPath{i});
        plot_blinks19(60, 150, 50, 100, s.CCDpixelx_SM, s.CCDpixely_SM, ...
            s.imageStacksPathByRuns{i}, s.calFileSMPath{i}, s.calFileSM_Z_Wxy{i}, ...
            s.cellDescrip, s.date, s.runNum(i), s.calFileSM_Rxy_Z{i}, n_image, s.vispPath, s.imagingDepthInMicron(i));  % 50, 150
        % xyresolution, zresolution, pixsize, zdist  (10, 40, 10, 40)

        %close all;
        t_dstorm = toc(t_dstorm)/60;
        cprintf('*blue', 'Run#%g finisihed successfully in %.1f min. \n\n', s.runNum(i), t_dstorm);
    end
    
    % writing the analysis summary => s object and command window
    saveSummary(s);
      
end

