function plot_blinks19(xyresolution, zresolution, pixsize, zdist, CCDpixelx, CCDpixely, directoryName, directoryName_cal, filePath_ZWxy, imagetitle, DateRec, RunNo, Rxy_Z_table, nimage, vispPath, depth, drawtype, frameGroup)
%%
% - Additional exclusion criteria
% - Two options for reading position data file, i.e. position.dat and position.mat
% - Copying the PreSTORM images into the ViSP files path

%global A Xc Yc Zc dX dY dZ wX wY Rxy dwX dwY R frameNum

for n_plot=2:2
    %
    %filePath_Rxy = fullfile(directoryName,'\Rxy.dat');
    ZWxy = filePath_ZWxy;
    %ZWxy = load(filePath_ZWxy); % 1:wx, 2:wy 3:Z
    %
    PlotSquare = false;
    PlotPoint = false;
    PlotError = true;
    PlotSingles = true;
    
    if nargin < 16
        imagetitle = 'STORM';
    end
    
    if nargin > 16
        wordflag = 0;
        if (numel(strfind(drawtype,'no_singles')) > 0)
            wordflag = wordflag + 1;
            PlotSingles = false;
        end
        if (numel(strfind(drawtype,'default')) > 0)
            wordflag = wordflag + 1;
            PlotSingles = true;
        end
        if (numel(strfind(drawtype,'plot_multiples')) > 0)
            wordflag = wordflag + 1;
            PlotSingles = false;
        end
        if (numel(strfind(drawtype,'plot_singles')) > 0)
            wordflag = wordflag + 1;
            PlotSingles = true;
        end
        if(numel(strfind(drawtype,'plot_square')) > 0)
            wordflag = wordflag + 1;
            PlotSquare = true;
            PlotError = false;
            PlotPoint = false;
        end
        if(numel(strfind(drawtype, 'plot_point')) > 0)
            wordflag = wordflag + 1;
            PlotSquare = false;
            PlotError = false;
            PlotPoint = true;
        end
        if(numel(strfind(drawtype, 'plot_error')) > 0)
            wordflag = wordflag + 1;
            PlotSquare = false;
            PlotError = true;
            PlotPoint = false;
        end
        nwords = length(regexp(drawtype, '\s+'))+1;
        if (nwords ~= wordflag)
            disp(num2str(nwords))
            disp(['Unknown option in drawing type: ' drawtype ]);
        end
    end
    
    if nargin < 17
        frameGroup = 20;
    end
    
    % Z limits based on the calibration curve
    %filePath_Rxy_Z_appRange = fullfile(directoryName_cal,'\Rxy_Z_appRange.dat');
    %Rxy_Z_appRange = load(filePath_Rxy_Z_appRange);
    llimZ = min(Rxy_Z_table(:,2));%min(Rxy_Z_appRange(:,2))-300;
    ulimZ = max(Rxy_Z_table(:,2));%max(Rxy_Z_appRange(:,2))+300;
    
    llimRxy = min(Rxy_Z_table(:,1))-1;%min(Rxy_Z_appRange(:,1))-0.1;
    ulimRxy = max(Rxy_Z_table(:,1))+1;% max(Rxy_Z_appRange(:,1));
    
    ulimWx = max(ZWxy(:,1));
    ulimWy = max(ZWxy(:,2));
    
    disp(['dSTORM - Analysing ' directoryName ]);
    disp([' Pixel size = ' num2str(pixsize) 'nm']);
    disp([' XY Resolution = ' num2str(xyresolution) 'nm']);
    disp([' Z Resolution = ' num2str(zresolution) 'nm']);
    disp([' Z Spacing = ' num2str(zdist) 'nm']);
    fprintf('\n');
    %
    % Define STORM pixel size
    xpix=pixsize; % nm %<=======================
    ypix=xpix; % nm %<=======================
    
    % define PSF in STORM (delta is the width
    %delta=7; %(nm) width of pixel %<=======================
    % Limits
    llimR = 0.3;
    ulimw = 3.2; %3
    llimw = 0.5; %1
    ulimA = 500;%150;%600;
    llimA = 2;%25;%25;  2
    llimdxy = 0.001;
    ulimdx = xyresolution/CCDpixelx;
    ulimdy = xyresolution/CCDpixely;
    llim_dz = 5;
    %     ulimdwx = zresolution/(sqrt(2)*CCDpixelx);
    %     ulimdwy = zresolution/(sqrt(2)*CCDpixely);
    %
    % 3D fitting Z polynomial coef.: Coef(1,:)=Rxy-Z; Coef(2,:)=Rxy-X; Coef(3,:)=Rxy-Y
    % Coef = load(CalFileName);
    % load position (pixels)
    cprintf('Reading the position data ... \n\n');
    filePath = fullfile(directoryName, 'position.mat');
    if exist(filePath,'file') == 2
        load(filePath)
    else
        filePath = fullfile(directoryName, 'position.dat');
        position = load(filePath);
    end
    %
    % plot wide-field image
    xmaxw = NaN;
    ymaxw = NaN;
    
    %filePath = fullfile(directoryName,'WideField.tif');
    %if exist(filePath,'file')
    %    info = imfinfo(filePath);
    %    xmaxw = info.Width ;
    %    ymaxw = info.Height;
    %else
    filePath = fullfile(directoryName,'WideField.tif');  %< input STORM images +++++++
    if exist(filePath,'file') == 2
        %info = imfinfo(filePath);
        xmaxp = 512; %info.Width;
        ymaxp = 512; %info.Height;
        %    WFI = imread(I, 1, 'Info', info);
        %    figure(1),subplot(211),imshow(WFI,[min(WFI(:)) max(WFI(:))]), title('Wide-Field Image');
    else
        error('Widefield pre-STORM image does not exist, its size is required for SM image recunstruction.  => WideFiled.tif \n');
    end
    %end
    xmax = min(xmaxw, xmaxp);
    ymax = min(ymaxw, ymaxp);
    % create arrays
    % fill arrays with position values
    A = position(:,2);
    dX = position(:,9);
    dY = position(:,11);
    wX = position(:,4);
    wY = position(:,6);
    dwX = position(:,10);
    dwY = position(:,12);
    R = position(:,13);
    FG = position(:, 14);
    fitstatus = position(:,7);
    Rxy = wX./wY;
    
    %% Finding the best fit for Wy vs. Wx
    bestFittingCurve_option = 2; % 1: Yes, 2: No
    [nd,~] = size(position);
    
    if bestFittingCurve_option == 1
        
        t_tic = tic;
        cprintf('*blue', 'Finding an OPTIMUM CALIBRATION CURVE by fitting the blinks ...\n');
        k_count = 0;
        if nd > 10000;
            n_pos = 10000;
        else
            n_pos = nd;
        end
        selBlinks = zeros(n_pos,3);
        while k_count < n_pos
            %j = randi([1+1e3 nd-1e3], 1,1); % a random number to extract a point form position file
            j = randi([1 nd], 1,1);
            if (wX(j,1) < ulimWx) && (wX(j,1) > llimw)
                if (wY(j,1) < ulimWy) && (wY(j,1) > llimw)
                    if (R(j) < 1) && (R(j) > llimR)
                        if (A(j) > llimA) && (A(j) < ulimA)
                            %if (dwX(j) < ulimdwx) && (dwY(j) < ulimdwy)
                            k_count = k_count + 1;
                            selBlinks(k_count,1) = wX(j,1);
                            selBlinks(k_count,2) = wY(j,1);
                            %end
                        end
                    end
                end
                
            end
        end
        selBlinks = selBlinks(any(selBlinks,2),:);
        [Wx_Wy_Z] = determineBestWyWxCurve(selBlinks, Rxy_Z_table, directoryName, ZWxy);
        Rxy_Z_table_bestFit = [Wx_Wy_Z(:,4) Wx_Wy_Z(:,3)];
        
        filePath_pos = [directoryName,'\Wx_Wy_Z.dat'];
        save(filePath_pos,'Wx_Wy_Z','-ASCII');
        
        h = figure;
        plot(ZWxy(:,1), ZWxy(:,2), '.-r', Wx_Wy_Z(:,1), Wx_Wy_Z(:,2), '.-b');
        legend('beads', 'blinks')
        xlabel('Wx'); ylabel('Wy');
        axis([0 5 0 5]); axis square
        
        fid = [directoryName '\calCurve_beadsVsBlinks.fig'];
        saveas(h, fid);
        
        t_toc = toc(t_tic);
        fprintf('Elapsed time for finding an optimum cal curve is %.1f min\n\n', t_toc/60);
        
    else
        fprintf('Disabled: the search for OPTIMUM CALIBRATION CURVE by fitting the blinks ...\n\n');
        calFile_path_bestCalCurveOff = fullfile(directoryName_cal,'ZWxy.dat');
        [Wx_Wy_Z] = load(calFile_path_bestCalCurveOff);
        filePath_pos = [directoryName,'\Wx_Wy_Z.dat'];
        save(filePath_pos,'Wx_Wy_Z','-ASCII');
    end
    
    %%
    figure_WxWy = figure;
    xedges = linspace(0,3,200); yedges = linspace(0,3,200);
    histmat = hist2(wX, wY, xedges, yedges);
    pcolor(xedges,yedges,histmat'); colorbar ; axis square tight;
    xlabel('wX');
    ylabel('wY');
    
    hold on;
    plot(Wx_Wy_Z(:,1), Wx_Wy_Z(:,2), '.k');
    axis([0 3 0 3])
    % ZWxy = load(filePath_ZWxy);
    % plot(ZWxy(:,1), ZWxy(:,2), '.k');
    filePath_WxWy = [directoryName,'\figure_WxWy.jpg'];
    saveas(figure_WxWy ,filePath_WxWy,'jpg');
    filePath_WxWy = [directoryName,'\figure_WxWy.fig'];
    saveas(figure_WxWy ,filePath_WxWy);
    
    %% calculate Z (nm) and shift x & y (pixel)
    %
    % Calculate Z using max(wx, wy)
    cprintf('*blue', 'Z-Localization ... \n');
    if bestFittingCurve_option == 1
        [Xc,Yc,Zc,dZ, distance]  = cal_xyz14(position(:,3),position(:,4),dwX(:,1), position(:,5),position(:,6),dwY(:,1), Wx_Wy_Z, Rxy_Z_table); % xc, Wx, dwX,yc,Wy,dwY,ZWxy
    else
        [Xc,Yc,Zc,dZ, distance]  = cal_xyz14(position(:,3),position(:,4),dwX(:,1), position(:,5),position(:,6),dwY(:,1), ZWxy, Rxy_Z_table); % xc, Wx, dwX,yc,Wy,dwY,ZWxy
    end
    
    %[Zc,dZ] = Z_loc_usingRxy02(wX, dwX, wY, dwY, Rxy_Z_table);
    
    fig_dist = figure;
    set(fig_dist,'name', 'Distance from cal curve');
    [count, bins] = hist(distance,200);
    bar(bins, count, 'hist');
    %axis([0 20 0 max(count)]);
    xlabel('Distance (px)', 'FontSize', 16);
    ylabel('Frequency', 'FontSize', 16);
    
    Zc(distance > 4) = NaN;
    dZ(distance > 4) = NaN;
    wX(isnan(Zc))=NaN;
    wY(isnan(Zc))=NaN;
    dwX(isnan(Zc))=NaN;
    dwY(isnan(Zc))=NaN;
    
    %
    %[Zc,dZ] = Z_loc_usingRxy02(wX, dwX, wY, dwY, Rxy_Z_table);
    %[Zc,dZ] = Z_loc_usingRxy02(wX, dwX, wY, dwY, Rxy_Z_table_bestFit);

    %%
    
    % Offline Drift Corrections
    drift_correction = 0; % yes:1, no: 0
    if drift_correction == 1
        if n_plot == 2
            Xc=Xc-position(:,17);
            Yc=Yc-position(:,16);
            Zc=Zc-position(:,18);
        end
    end
    %
    %
    nx_pixel=xmax;
    ny_pixel=ymax;
    %% criteria for exclusion
    crit_dis = sum(isnan(Zc));
    crit_total = sum(isnan(Zc));
    
    Zc(dX > ulimdx) = NaN; % 
    crit_xres = sum(isnan(Zc)) - crit_total;
    crit_total = sum(isnan(Zc));
    
    Zc(dY > ulimdy) = NaN; % 
    crit_yres = sum(isnan(Zc)) - crit_total;
    crit_total = sum(isnan(Zc));
    
    %Zc(dZ > zresolution | dZ < 5) = NaN; % 100
    Zc(dZ > zresolution | dZ < llim_dz) = NaN;
    crit_zres = sum(isnan(Zc)) - crit_total;
    crit_total = sum(isnan(Zc));
    
    Zc(A < llimA | A > ulimA) = NaN;
    crit_amp = sum(isnan(Zc)) - crit_total;
    crit_total = sum(isnan(Zc));
    
    Zc(R > 1 | R < llimR ) = NaN;
    crit_R = sum(isnan(Zc)) - crit_total;
    crit_total = sum(isnan(Zc));
    
    %Zc(wX < llimw  | wX > ulimw ) = NaN;
    Zc(wX < llimw  | wX > ulimWx ) = NaN;
    crit_Wx = sum(isnan(Zc)) - crit_total;
    crit_total = sum(isnan(Zc));
    
    %Zc(wY < llimw | wY > ulimw ) = NaN;
    Zc(wY < llimw | wY > ulimWy ) = NaN;
    crit_Wy = sum(isnan(Zc)) - crit_total;
    crit_total = sum(isnan(Zc));
    
    Zc(Zc < llimZ | Zc > ulimZ ) = NaN;
    crit_Zlim = sum(isnan(Zc)) - crit_total;
    crit_total = sum(isnan(Zc));
    
    Zc(Rxy > ulimRxy | Rxy < llimRxy) = NaN;
    crit_Rxy = sum(isnan(Zc)) - crit_total;
    crit_total = sum(isnan(Zc));
    
    %     Zc(dwX < ulimdwx) = NaN;
    %     crit_dwx = sum(isnan(Zc)) - crit_total;
    %     crit_total = sum(isnan(Zc));
    %
    %     Zc(dwY < ulimdwy) = NaN;
    %     crit_dwy = sum(isnan(Zc)) - crit_total;
    %     crit_total = sum(isnan(Zc));
    
    
    A(isnan(Zc))=NaN;
    R(isnan(Zc))=NaN;
    wX(isnan(Zc))=NaN;
    wY(isnan(Zc))=NaN;
    dwX(isnan(Zc))=NaN;
    dwY(isnan(Zc))=NaN;
    Xc(isnan(Zc))=NaN;
    Yc(isnan(Zc))=NaN;
    dX(isnan(Zc))=NaN;
    dY(isnan(Zc))=NaN;
    dZ(isnan(Zc))= NaN;
    Rxy(isnan(Zc))= NaN;
    
    
    %% 
    fprintf('\n');
    cprintf('*blue', 'Summary\n\n');
    cprintf('*black','%d failed Z-distance-to-calibration-curve criteria. \n', crit_dis);
    cprintf('*black', '%d failed X resolution criteria. \n', crit_xres );
    cprintf('*black', '%d failed Y resolution criteria. \n', crit_yres);
    cprintf('*black','%d failed amplitude critera. \n',crit_amp);
    cprintf('*black','%d failed goodness of fit critera. \n', crit_R);
    cprintf('*black','%d failed Z resolution criteria. \n', crit_zres);
    cprintf('*black','%d failed Z range criteria. \n', crit_Zlim);
    cprintf('*black','%d failed Wx criteria. \n', crit_Wx);
    cprintf('*black','%d failed Wy criteria. \n', crit_Wy);
    cprintf('*black','%d failed Rxy criteria. \n', crit_Rxy);
    %     fprintf('%d failed dwx criteria\n', crit_dwx);
    %     fprintf('%d failed dwy criteria\n', crit_dwy);
    cprintf('*blue','total failures = %d  => %g%% of total blinks. \n\n', crit_total, round(100*crit_total/nd));
    cprintf('*black', 'Total blinks (started with): %d.\n', nd);
    %%
    % plot histograms
    figure_hist = figure;
    set(figure_hist,'name', 'Histograms');
    
    subplot(621);
    [count, bins] = hist(A,50);
    bar(bins, count, 'hist');
    axis([0 ulimA+100 0 max(count)]);
    %axis square;
    xlabel('A', 'FontSize', 16);
    
    subplot(622);
    [count, bins] = hist(dX*CCDpixelx,100);
    bar(bins, count, 'hist');
    axis([0 2*xyresolution 0 max(count)]);
    %axis square;
    xlabel('dX (nm)', 'FontSize', 16);
    
    subplot(623)
    [count, bins] = hist(dY*CCDpixely,100);
    bar(bins, count, 'hist');
    axis([0 2*xyresolution 0 max(count)]);
    %axis square;
    %axis([0 20 0 max(count)]);
    xlabel('dY (nm)', 'FontSize', 16);
    
    subplot(624)
    [count, bins] = hist(dZ,100);
    bar(bins, count, 'hist');
    axis([0 2*zresolution 0 max(count)]);
    %axis square;
    xlabel('dZ (nm)', 'FontSize', 16);
    
    subplot(625);
    [count, bins] = hist(wX,100) ;
    bar(bins, count, 'hist');
    axis([llimw ulimw 0 max(count)]);
    %axis square;
    xlabel('Wx (px)', 'FontSize', 16);
    
    subplot(626);
    [count, bins] = hist(dwX,100) ;
    bar(bins, count, 'hist');
    axis([0 0.5 0 max(count)]);
    %axis square;
    xlabel('dWx (px)', 'FontSize', 16);
    
    subplot(627);
    [count, bins] = hist(wY,100) ;
    bar(bins, count, 'hist');
    axis([llimw ulimw 0 max(count)]);
    %axis square;
    xlabel('Wy (px)', 'FontSize', 16);
    
    subplot(628);
    [count, bins] = hist(dwY,100) ;
    bar(bins, count, 'hist');
    axis([0 0.5 0 max(count)]);
    %axis square;
    xlabel('dWy (px)', 'FontSize', 16);
    
    subplot(629);
    [count, bins] = hist(Rxy,100) ;
    bar(bins, count, 'hist');
    axis([0 3 0 max(count)]);
    %axis square;
    xlabel('Rxy', 'FontSize', 16);
    
    subplot(6,2,10);
    [count, bins] = hist(Zc,100) ;
    bar(bins, count, 'hist');
    axis([-1000 1000  0 max(count)]);
    %axis square;
    xlabel('Zc (nm)', 'FontSize', 16);
    
    subplot(6,2,11);
    [count, bins] = hist(R,100) ;
    bar(bins, count, 'hist');
    axis([llimR/2 1  0 max(count)]);
    %axis square;
    xlabel('Goodness of Fit', 'FontSize', 16);
    
    %subplot(6,2,12);
    
    
    filePath_hist = [directoryName,'\Plot_Hist.jpg'];
    saveas(figure_hist ,filePath_hist,'jpg');
    
    %%
    figure_Z_distance = figure;
    set(figure_Z_distance, 'name', 'Distance from cal curve_histogram');
    xedges = linspace(-1000,1000,100); yedges = linspace(0,0.3,100);
    histmaty = hist2(Zc, distance, xedges, yedges);
    pcolor(xedges,yedges,histmaty'); colorbar ; axis square tight;
    xlabel('z (nm)', 'FontSize', 16);
    ylabel('dsitance (pixel)', 'FontSize', 16);
    filePath_Z_distance = [directoryName,'\Plot_Z_distance.jpg'];
    saveas(figure_Z_distance ,filePath_Z_distance,'jpg')
    
    %
    figure_Z_wxy = figure;
    set(figure_Z_wxy, 'name', 'WxWy - Z');
    xedges = linspace(-1000,1000,100); yedges = linspace(0,4,100);
    histmatx = hist2(Zc, wX, xedges, yedges);
    % pcolor(xedges,yedges,histmatx'); colorbar ; axis square tight;
    % xlabel('Z');
    % ylabel('wX');
    % figure(54)
    histmaty = hist2(Zc, wY, xedges, yedges);
    % pcolor(xedges,yedges,histmaty'); colorbar ; axis square tight;
    % xlabel('Z');
    % ylabel('wY');
    histmatxy=histmatx+histmaty;
    pcolor(xedges,yedges,histmatxy'); colorbar ; axis square tight;
    xlabel('Z');
    ylabel('wX/wY');
    filePath_Z_wxy = [directoryName,'\Plot_Z_wxy.jpg'];
    saveas(figure_Z_wxy ,filePath_Z_wxy,'jpg')
    %% Plotting Z vs dX. dY
    
    figure_Z_dX_dY = figure;
    set(figure_Z_dX_dY, 'name', 'dX_dY vs. Z _ Hist')
    subplot(1,2,1);
    xedges = linspace(-1000,1000,100); yedges = linspace(0,20,100);
    histmaty = hist2(Zc, dX*CCDpixelx, xedges, yedges);
    pcolor(xedges,yedges,histmaty'); colorbar ; axis square tight;
    xlabel('z (nm)','FontSize', 16);
    ylabel('dXc (nm)','FontSize', 16);
    
    subplot(1,2,2);
    histmaty = hist2(Zc, dY*CCDpixely, xedges, yedges);
    pcolor(xedges,yedges,histmaty'); colorbar ; axis square tight;
    xlabel('z (nm)','FontSize', 16);
    ylabel('dYc (nm)','FontSize', 16);
    filePath_Z_dX_dY = [directoryName,'\Plot_Z_dX_dY.jpg'];
    saveas(figure_Z_dX_dY ,filePath_Z_dX_dY,'jpg');
    
    %
    %%
        %%
    figure_dWxdWydZ = figure;
    set(figure_dWxdWydZ, 'name', 'WxWy_Z_Hist');
    subplot(221);
    xedges = linspace(0,4,100); yedges = linspace(0,4,100);
    histmat = hist2(wX, wY, xedges, yedges);
    pcolor(xedges,yedges,histmat'); colorbar ; axis square tight;
    xlabel('wX','FontSize',16);
    ylabel('wY','FontSize',16);
    
    hold on;
    plot(Wx_Wy_Z(:,1), Wx_Wy_Z(:,2), '.k');
    axis([0 4 0 4]);
    
    subplot(222)
    xedges = linspace(1,3,100); yedges = linspace(0,0.5,100);
    histmat = hist2(wX, dwX, xedges, yedges);
    save 'histmat.mat' histmat
    pcolor(xedges,yedges,histmat'); colorbar ; axis square tight;
    xlabel('wX','FontSize',16);
    ylabel('dwX','FontSize',16);
    subplot(223)
    xedges = linspace(1,3,100); yedges = linspace(0,0.5,100);
    histmat = hist2(wY, dwY, xedges, yedges);
    pcolor(xedges,yedges,histmat'); colorbar ; axis square tight;
    xlabel('wY','FontSize',16);
    ylabel('dwY','FontSize',16);
    subplot(224)
    xedges = linspace(-1000,1000,100); yedges = linspace(0,150,100);
    histmat = hist2(Zc,dZ, xedges, yedges);
    pcolor(xedges,yedges,histmat'); colorbar ; axis square tight;
    xlabel('Zc','FontSize',16);
    ylabel('dZ','FontSize',16);
    filePath_WxWy = [directoryName,'\figure_dWxdWydZ.jpg'];
    saveas(figure_dWxdWydZ ,filePath_WxWy,'jpg');
    filePath_WxWy = [directoryName,'\figure_dWxdWydZ.fig'];
    saveas(figure_dWxdWydZ ,filePath_WxWy);
    
    %% Output VISP File
    Visp_output = 1; % 1: yes, 0: no
    
    if Visp_output == 1
        
        if n_plot == 1
            vispfile = [ directoryName '\' imagetitle '_' DateRec '_' num2str(RunNo) '_' num2str(depth) 'um' '.3dlp'];
        end
        if n_plot == 2
            vispfile = [ directoryName '\' imagetitle '_' DateRec '_' num2str(RunNo) '_' num2str(depth) 'um' '.3dlp'];
        end
        fprintf('Writing visp file %s ...',vispfile);
        fid_3dlp = fopen(vispfile,'w');
        % translate x-y pixels to nm
        Xnm=Xc.*CCDpixelx;
        Ynm=Yc.*CCDpixely;
        dXnm=dX.*CCDpixelx;
        dYnm=dY.*CCDpixely;
%         
%         Xnm = Xc;%.*CCDpixelx;
%         Ynm = Yc;%.*CCDpixely;
%         dXnm = dX;%.*CCDpixelx;
%         dYnm = dY;%.*CCDpixely;
        
        m3d=0;
        
        for i = 1:nd
            Nofile=position(i,14)-1;
            NumFrame= sum(nimage(1:Nofile))+position(i,15);
            data_visp(1,1) = Xnm(i,1);
            data_visp(1,2) = Ynm(i,1);
            data_visp(1,3) = Zc(i,1);
            data_visp(1,4) = dXnm(i,1);
            data_visp(1,5) = dYnm(i,1);
            data_visp(1,6) = dZ(i,1);
            data_visp(1,7) = A(i,1);
            flag_visp = sum(isnan(data_visp));
            if flag_visp == 0
                m3d=m3d+1;
                fprintf(fid_3dlp,'%7.3f \t %7.3f \t %7.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f \t %d \n',...
                    Xnm(i,1), Ynm(i,1), Zc(i,1), dXnm(i,1), dYnm(i,1), dZ(i,1), A(i,1), NumFrame);
                %     fprintf('%7.0f \t %7.0f \t %7.0f \t %6.1f \t %6.1f \t %6.1f \t %6.2f \t %d \n',...
                %         Xnm(i,1), Ynm(i,1), Zc(i,1), dXnm(i,1), dYnm(i,1), dZ(i,1), A(i,1), NumFrame);
            end
        end
        fprintf('%g points written \n',m3d);
        fclose(fid_3dlp);
        
        %    Visp_path = ['D:\' DateRec '\VISP' '\' num2str(RunNo)];
        %     if exist(Visp_path,'dir')~=7
        %         parentFolder =['D:\' DateRec '\VISP'];
        %         folderName=num2str(RunNo);
        %         mkdir(parentFolder,folderName);
        %     end
        vispfileCopy = [vispPath '\' imagetitle '_' DateRec '_' num2str(RunNo) '_' num2str(depth) 'um' '.3dlp'];
        copyfile(vispfile, vispfileCopy);
        
        PreSTORM_path = [directoryName '\PreSTORM.tif'];
        PreSTORM_copy = [vispPath '\PreSTORM_' num2str(RunNo) '.tif'];
        copyfile(PreSTORM_path, PreSTORM_copy);
        
    end

    
    %%
    %
    %Zc = Zc -min(Zc); % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % z range
    nz1 = llimZ; % nm
    %nz2=ceil((ulimZ-min(Zc))/100)*100; % nm
    nz2 = ulimZ; % nm
    % Calculate pixel ratio (regular/STORM)
    Rpx=CCDpixelx/xpix;
    Rpy=CCDpixely/ypix;
    % calculate number of pixel in STORM images
    nx = ceil(ny_pixel*Rpx);
    ny = ceil(nx_pixel*Rpy);
    nz = ceil((nz2-nz1)/zdist);
    
    Zc = Zc - llimZ;
    
    disp('Creating STORM image structures');
    disp(['Size :' num2str(ny) 'x' num2str(nx) 'x' num2str(nz)]);
    
    % Create STORM images
    image=zeros(nx,ny,nz,'single');
    %
    % filter
    disp('Identifying and filling in blink positions ...');
    nplot=0;
    acrit = 0;
    rcrit = 0;
    xyrescrit = 0;
    zrescrit = 0;
    pcrit = 0;
    pxycrit = 0;
    scrit = 0;
    outgroup = 0;
    fiterr = 0;
    %disp(['Requesting ' num2str(nd) ' int16 ']);
    Xp=zeros(nd, 1);
    Yp=zeros(nd, 1);
    Zp=zeros(nd, 1);
    A0=200; % amplitude (arbitrary but make sure the image will not saturate)
    for k=1:nd
        %criteria for exclusion
        %      if isnan(fitstatus(k,1))
        %         fiterr = fiterr + 1;
        %         continue;
        %      end
        %      if FG(k,1) > frameGroup
        %         outgroup = outgroup + 1;
        %         continue;
        %      end
        %      if(R(k,1) < llimR)
        %         rcrit = rcrit +1;
        %         continue;
        %      end
        %      if(wX(k,1) < llimw || wX(k,1) > ulimw || wY(k,1) < llimw || wY(k,1) > ulimw)
        %         scrit = scrit + 1;
        %         continue;
        %      end
        %      if(A(k,1) < llimA) || (A(k,1) > ulimA)
        %         acrit = acrit + 1;
        %         continue;
        %      end
        %      if isnan(Zc(k,1))
        %         pcrit = pcrit + 1;
        %         continue;
        %      end
        %      if( Xc(k,1) <= 0 || Yc(k,1) <= 0)
        %         pxycrit = pxycrit + 1;
        %         continue;
        %      end
        %      if (dX(k,1)  < llimdxy || dX(k,1) > ulimdx || dY(k,1) < llimdxy || dY(k,1) > ulimdy)
        %         xyrescrit = xyrescrit + 1;
        %         continue;
        %      end
        %      if dZ(k,1) > zresolution
        %         zrescrit = zrescrit + 1;
        %         continue;
        %      end
        % convert Z (nm) to pixel
        NanBool = isnan(Xc(k)) + isnan(Yc(k)) + isnan(Zc(k));
        if NanBool == 0
            Zpos=max(1,round(Zc(k,1)/zdist));
            % translate x-y to STORM pixels
            Xpos=round(Xc(k,1)*Rpx); %center  (pixel)
            Ypos=round(Yc(k,1)*Rpy); %center
            if(Xpos <= nx  && Ypos <= ny)
                image(Xpos, Ypos, Zpos) = image(Xpos,Ypos,Zpos)+A0;
                nplot=nplot+1;
                Xp(nplot, 1) = Xpos;
                Yp(nplot, 1) = Ypos;
                Zp(nplot, 1) = Zpos;
            end
        end
    end
    
    
    disp('Writing position_table ...');
    frameNum = zeros(nd,1);
    for i = 1:nd
        Nofile=position(i,14)-1;
        frameNum(i,1)= sum(nimage(1:Nofile))+position(i,15);
    end
    frameNum(isnan(Zc),1) = NaN;
    
    %     position_table = table(A, Xc, Yc, Zc, dX, dY, dZ, wX, wY, Rxy, dwX, dwY, R, frameNum);
    %     position_table.Properties.VariableDescriptions{'R'} = 'goodnessOfFit';
    %     position_table(isnan(Zc),:) = {NaN};
    %     %position_table = position_table(all(~isnan(position_table),2),:);
    %     filePath = fullfile(directoryName, 'position_table.csv');
    %     writetable(position_table, filePath);
    
    position_table = [A Xc Yc Zc+llimZ dX dY dZ wX wY Rxy dwX dwY R frameNum];
    position_table(isnan(Zc),:) = NaN;
    position_table = position_table(all(~isnan(position_table),2),:);
    filePath = fullfile(directoryName, 'position_table.mat');
    save( filePath, 'position_table');
    
    
    position(isnan(Zc),:) = NaN;
    position = position(all(~isnan(position),2),:);
    position_filtered = position;
    clear position;
    filePath = fullfile(directoryName, 'position_filtered.mat');
    save(filePath, 'position_filtered');
    
    
    clear position_filtered;
    clear A;
    clear R;
    clear dX;
    clear dY;
    clear dZ;
    clear dwX;
    clear dwY;
    clear wX;
    clear wY;
    
    
    %fprintf('%d failed shape criteria\n', scrit);
    %fprintf('%d failed XY position criteria\n', pxycrit);
    %fprintf('%d had fit errors\n',fiterr);
    %fprintf('%d excluded because frame group > %d\n',outgroup, frameGroup);
    
    imagetype = '';
    description = '';
    % if ~(PlotSingles)
    %   for m=1:nplot
    %     i = Xp(m, 1);
    %     j = Yp(m, 1);
    %     k = Zp(m, 1);
    %     if(image(i,j,k) == A0)
    %        image(i, j, k) = 0;
    %     end
    %   end
    %   imagetype = '_nosingle';
    %   description = 'No single blinks\n';
    % end
    % if PlotError
    %   imagetype = [imagetype  '_gauss_'];
    %   description = [description 'Plot is Gaussian Error\n'];
    % elseif PlotSquare
    %   imagetype = [imagetype '_square_'];
    %   description = [description 'Plot is 30 nm square\n'];
    % else
    %   imagetype = [imagetype '_point_'];
    %   description = [description 'Plot is Gaussian Error\n'];
    % end
    %
    % if ~PlotPoint
    %   sX=((30/xpix)-1)/2;
    %   sY=((30/ypix)-1)/2;
    %   wX0=delta/xpix;%dX(k,1)*Rpx;  %pixel
    %   wY0=delta/ypix;%dY(k,1)*Rpy*1; %pixel
    % %     wZ0=delta;%dZ0; % nm
    %   for m=1:nplot
    %     i = Xp(m,1);
    %     j = Yp(m,1);
    %     k = Zp(m,1);
    %     if(image(i,j,k) > 0)
    %         A1 =image(i,j,k);
    %         if(PlotError)
    %    % calculate the area to place the PSF
    %            x1=floor((i-4*wX0));
    %            x2=ceil((i+4*wX0));
    %            y1=floor((j-4*wY0));
    %            y2=ceil((j+4*wY0));
    %         else
    %            x1=floor(i-sX);
    %            x2=ceil(i+sX);
    %            y1=floor(j-sY);
    %            y2=ceil(j+sY);
    %         end
    %       % avoid edges
    %         if x1 < 1
    %            x1=1;
    %         end
    %         if x2 > nx
    %            x2=nx;
    %         end
    %         if y1 < 1
    %            y1=1;
    %         end
    %         if y2 > ny
    %            y2=ny;
    %         end
    %       % Place the PSF
    %         kz = k;
    %         if PlotError
    %            for ky=y1:y2
    %               yexp = (ky-j)/wY0;
    %               FY=exp(-(yexp*yexp));
    %               for kx=x1:x2
    %                  xexp = (kx-i)/wX0;
    %                  if((xexp + yexp) ~= 0)
    %                    image(kx,ky,kz)= image(kx,ky,kz)+A1*exp(-(xexp*xexp))*FY;
    %                  end
    %               end
    %            end
    %         end
    %         if PlotSquare
    %            for ky=y1:y2
    %               for kx=x1:x2
    %                  if (kx ~= i && ky ~= j)
    %                    image(kx,ky,kz)= image(kx,ky,kz)+A1;
    %                  end;
    %               end
    %            end
    %         end %..if PlotError
    %      end % if image > 0
    %   end % for l
    % end
    
    %
    %
    addtitle ='';
    if frameGroup < 20
        addtitle = ['fg' num2str(frameGroup) '_'];
    end
    if n_plot == 1
        fn3d = [imagetitle addtitle imagetype 'pix_' num2str(xpix) 'nm_z_' num2str(zdist) 'nm_xyres_' num2str(xyresolution) 'nm_zres_' num2str(zresolution) 'nm_R.tif'];
        fprintf('Writing image %s ...',fn3d);
    end
    % if n_plot == 2
    % fn3d = [imagetitle addtitle imagetype 'pix_' num2str(xpix) 'nm_z_' num2str(zdist) 'nm_xyres_' num2str(xyresolution) 'nm_zres_' num2str(zresolution) 'nm_W.tif'];
    % fprintf('Writing image %s ...',fn3d);
    % end
    
    if n_plot == 2
        fn3d = [imagetitle addtitle '_gauss_' 'pix_' num2str(xpix) 'nm_z_' num2str(zdist) 'nm_xyres_' num2str(xyresolution) 'nm_zres_' num2str(zresolution) 'nm.tif'];
        fprintf('Writing image %s ...',fn3d);
    end
    
    s=Tiff(fullfile(directoryName,fn3d),'w');
    iminf.Software='dSTORM UBC';
    %iminf.ImageDescription = [ 'dSTORM\nPixel size = ' num2str(xpix) 'nm\n' description '\nXY Resolution = ' num2str(xyresolution) 'nm\nZ Resolution = ' num2str(zresolution) 'nm\nZ Spacing = ' num2str(zdist)];
    iminf.ImageLength=nx;
    iminf.ImageWidth=ny;
    %iminf.ResolutionUnit=Tiff.ResolutionUnit.Centimeter;
    %iminf.XResolution=1e7/xpix;
    %iminf.YResolution=iminf.XResolution;
    iminf.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    iminf.DateTime=datestr(now,'yyyy:mm:dd HH:MM:SS');
    iminf.SamplesPerPixel=1;
    iminf.BitsPerSample=16;
    iminf.Photometric=Tiff.Photometric.MinIsBlack;
    for ks=1:nz
        s.setTag(iminf);
        s.write(uint16(image(:,:,ks)));
        s.writeDirectory();
        % fprintf('.');
    end
    fprintf('\n');
    s.close();
    clear s;
    clear image;
    
    cprintf('*red','Blinks plotted : %d.\n\n',nplot);
    cprintf('*String', 'Z-localization DONE! \n\n');
    
end