function [Drift, nframe, bead_xyz]  = drift_bead(directoryName, nbead, Coef,CCDpixely, CCDpixelx, CCDpixely_newton, CCDpixelx_newton, nimage)
%global Drift
%************************************
if nbead == 0
    nframe=sum(nimage(:));
    Drift=ones(nframe,5);
    bead_xyz=0; %(1:x, 2:y, 3:z, 4:dx, 5:dy, 6:dz)
elseif nbead > 0
    %
    % Load data files
    filePath = fullfile(directoryName,'drift.mat');
    Volt=cell2mat(struct2cell(load(filePath)));
    filePath = fullfile(directoryName,'newtonTime.mat');
    NewtonTime=cell2mat(struct2cell(load(filePath)));
    %
    % Set new time zero
    Tn(:,1)=NewtonTime(:,3)*24*60*60+NewtonTime(:,4)*60*60+NewtonTime(:,5)*60+NewtonTime(:,6);
    Tn(:,2)=NewtonTime(:,10)*24*60*60+NewtonTime(:,11)*60*60+NewtonTime(:,12)*60+NewtonTime(:,13);
    filePath = fullfile(directoryName,'iXonTime.mat');
    iXonTime=cell2mat(struct2cell(load(filePath)));
    Tx(:,1)=iXonTime(:,4)*24*60*60+iXonTime(:,5)*60*60+iXonTime(:,6)*60+iXonTime(:,7);
    Tx=Tx-Tn(1,1); % reset time=0
    Tn=Tn-Tn(1,1); % reset time=0
    %
    n1=sum(nimage(:));% find number of frames from iXon
    [n2,n3]=size(Tx);
    nframe=min(n1,n2);
    Drift=zeros(nframe,5);% (flag, Frame_newton, x, y, z)
    %
    filePath_bead=[directoryName,'\Bead',int2str(1),'.tif'];
    I = (filePath_bead);
    info = imfinfo(I);
    Nimage_bead = numel(info);
    position_bead=zeros(Nimage_bead,14,nbead);
    bead_xyz=zeros(Nimage_bead,6,nbead); %(1:x, 2:y, 3:z, 4:dx, 5:dy, 6:dz)
    frame_newton=1;
    %
    % identify newton frame number for iXon
    for frame_ixon=1:nframe
      if Tx(frame_ixon,1) >= 0
        while Tx(frame_ixon,1) > Tn(frame_newton,1)
            frame_newton=frame_newton+1;
        end
        frame_newton=frame_newton-1;
        if Tx(frame_ixon,1) >= Tn(frame_newton,1) && Tx(frame_ixon,1) < Tn(frame_newton,2)
            Drift(frame_ixon,1)=1; %0
            Drift(frame_ixon,2)=frame_newton;
        elseif Tx(frame_ixon,1) > Tn(frame_newton,1) && Tx(frame_ixon,1) >= Tn(frame_newton,2)
            Drift(frame_ixon,1)=1;
            Drift(frame_ixon,2)=frame_newton;
        else
            Drift(frame_ixon,1)=1; %0
            Drift(frame_ixon,2)=frame_newton;
        end
      else
        Drift(frame_ixon,1)=1; %0
        Drift(frame_ixon,2)=frame_newton;
      end
%         fprintf('frame_ixon: %g; Tx: %g, Tn %g - %g \n',frame_ixon, Tx(frame_ixon,1),Tn(frame_newton,1:2));
%         fprintf('Drift: %g, %g \n',Drift(frame_ixon,1),Drift(frame_ixon,2));
%          pause(1e-2)
    end
    %
    % Import bead images and fitting
    for n=1:nbead
       filePath_bead=[directoryName,'\Bead',int2str(n),'.tif'];
       I = (filePath_bead);
       info = imfinfo(I);
       Nimage_bead = numel(info);
       temp = imread(I, 1, 'Info', info);
       [nx_bead, ny_bead]=size(temp);
       imagbead=zeros(nx_bead, ny_bead, Nimage_bead);
       for k=1:Nimage_bead
           imagbead(:,:,k) = transpose(imread(I, k, 'Info', info));
       end
       fprintf('Fitting bead No. %g ... ',n);
       %
       [P]  = fit_erf_only(imagbead);
       position_bead(:,:,n)=P; 
       %[1:6]:a0, a, x0, wx, y0, wy
       %[7:12]:da0, da, dx0, dwx, dy0, dwy
       %[13,14]: R, k
    end
    %
    % Save position data
    filePath_position_bead = [directoryName,'\position_bead.mat'];
    fid_position = fopen(filePath_position_bead,'w');
    save(filePath_position_bead,'position_bead','-mat')
    fclose('all');
    %
    % ploting
%     figure(10)
%     subplot(711), plot(1:Nimage_bead, position_bead(:,1,1),'.'),title('bg');
%     subplot(712), plot(1:Nimage_bead, position_bead(:,2,1),'.'),title('A');
%     subplot(713), plot(1:Nimage_bead, position_bead(:,3,1),'.'),title('x0');
%     subplot(714), plot(1:Nimage_bead, position_bead(:,4,1),'.'),title('wx');
%     subplot(715), plot(1:Nimage_bead, position_bead(:,5,1),'.'),title('y0');
%     subplot(716), plot(1:Nimage_bead, position_bead(:,6,1),'.'),title('wy');
%     subplot(717), plot(1:Nimage_bead, position_bead(:,13,1),'.'),title('R');
%     figure(11)
%     subplot(611), plot(1:Nimage_bead, position_bead(:,7,1),'.'),title('dbg');
%     subplot(612), plot(1:Nimage_bead, position_bead(:,8,1),'.'),title('dA');
%     subplot(613), plot(1:Nimage_bead, position_bead(:,9,1),'.'),title('dx0');
%     subplot(614), plot(1:Nimage_bead, position_bead(:,10,1),'.'),title('dwx');
%     subplot(615), plot(1:Nimage_bead, position_bead(:,11,1),'.'),title('dy0');
%     subplot(616), plot(1:Nimage_bead, position_bead(:,12,1),'.'),title('dwy'); 
    %
    % calculate Drift corection for dx, dy, dz
    for n=1:nbead
        [Xc,Yc,Zc]  = cal_xyz03(position_bead(:,3,n),position_bead(:,4,n),...
                            position_bead(:,5,n),position_bead(:,6,n),Coef);
        bead_xyz(:,1,n)=Xc; % x0 (pixel)
        bead_xyz(:,2,n)=Yc; % y0 (pixel)
        bead_xyz(:,3,n)=Zc; % z0 (nm)
        bead_xyz(:,4,n)=position_bead(:,9,n); % dx0  (pixel)
        bead_xyz(:,5,n)=position_bead(:,11,n); % dy0 (pixel)
        % Error Propagation: dRxy/Rxy=sqrt((dwX/wX)^2+dwY/wY)^2)
        wX=position_bead(:,4,n); 
        dwX=position_bead(:,10,n);
        wY=position_bead(:,6,n);
        dwY=position_bead(:,12,n);
        Rxy=wX./wY;
        dRxy =sqrt((dwX./wX).^2+(dwY./wY).^2).*Rxy;
        %Calculate slope (Z vs. Rxy)
        Slope_Z0 =Coef(1,2)+Rxy.*(2*Coef(1,3)+Rxy.*(3*Coef(1,4)+Rxy.*(4*Coef(1,5)+Rxy.*(5*Coef(1,6)+Rxy*6*Coef(1,7)))));  %slope
        % Convert dRxy to dZ
        bead_xyz(:,6,n)=abs(Slope_Z0.*dRxy); % dz0(nm)
    end
    % remove bad fittings
    for n=1:nbead
        for k=2:Nimage_bead
            if sum(isnan(bead_xyz(k,:,n))) ~=0
                bead_xyz(k,1:6,n)=NaN;
            end
        end
        % ploting
%         figure
%         subplot(611), plot(1:Nimage_bead, bead_xyz(:,1,n),'.'),title('x0');
%         subplot(612), plot(1:Nimage_bead, bead_xyz(:,2,n),'.'),title('y0');
%         subplot(613), plot(1:Nimage_bead, bead_xyz(:,3,n),'.'),title('z0');
%         subplot(614), plot(1:Nimage_bead, bead_xyz(:,4,n),'.'),title('dx0');
%         subplot(615), plot(1:Nimage_bead, bead_xyz(:,5,n),'.'),title('dy0');
%         subplot(616), plot(1:Nimage_bead, bead_xyz(:,6,n),'.'),title('dz0'); 
    end
    % reset zeros
    flag0=0;
    for k=1:Nimage_bead
        if flag0 ==0
            temp01=isnan(bead_xyz(k,:,:));
            if sum(temp01(:)) == 0
               flag0=k;
            end
        end
    end
    %
    for n=1:nbead
        bead_xyz(:,1,n)=bead_xyz(:,1,n)-bead_xyz(flag0,1,n);
        bead_xyz(:,2,n)=bead_xyz(:,2,n)-bead_xyz(flag0,2,n);
        bead_xyz(:,3,n)=bead_xyz(:,3,n)-bead_xyz(flag0,3,n);
    end
    % save
    filePath_position_bead = [directoryName,'\bead_xyz.mat'];
    fid_position = fopen(filePath_position_bead,'w');
    save(filePath_position_bead,'bead_xyz','-mat')
    fclose('all');
%     for n=1:nbead 
%         figure
%         subplot(611), plot(1:Nimage_bead, bead_xyz(:,1,n),'.'),title('x0');
%         subplot(612), plot(1:Nimage_bead, bead_xyz(:,2,n),'.'),title('y0');
%         subplot(613), plot(1:Nimage_bead, bead_xyz(:,3,n),'.'),title('z0');
%         subplot(614), plot(1:Nimage_bead, bead_xyz(:,4,n),'.'),title('dx0');
%         subplot(615), plot(1:Nimage_bead, bead_xyz(:,5,n),'.'),title('dy0');
%         subplot(616), plot(1:Nimage_bead, bead_xyz(:,6,n),'.'),title('dz0');
%     end
%     %
    % calculate final drift corrections for x, y, and z
    % Estimate the likelihood positions using weighting factor ~ inverse square of the error
    for n=1:nframe
        FN=Drift(n,2);
        sumD=0;
        partition=0;
        flag_Drift=0;
        for m=3:5
            sumD=0;
            partition=0;
            for k=1:nbead
                if (bead_xyz(FN,m+1,k) ~=0) && sum(isnan(bead_xyz(FN,:,k)))==0
                    sumD=sumD+bead_xyz(FN,m-2,k)/bead_xyz(FN,m+1,k)^2;
                    partition=partition+1/bead_xyz(FN,m+1,k)^2;
                end
            end
            if partition~=0
                Drift(n,m)=sumD/partition;% (1:flag, 2:Frame_newton, 3:x, 4:y, 5:z)
                flag_Drift=flag_Drift+1;
            end
        end
        if flag_Drift == 0 
            Drift(n,1)=1; %0
        end
    end
%     % Find 1st nonzeros frame
%     n_nz=1;
%     for n=1:nframe
%         if Drift(n_nz,1) ~= 0
%            n_nz=n;
%            break
%         end 
%     end
%     % Reset zeros
%     Drift(n_nz:nframe,3)=(Drift(n_nz:nframe,3)-Drift(n_nz,3))*CCDpixely_newton/CCDpixelx;% (1:flag, 2:Frame_newton, 3:x, 4:y, 5:z)
%     Drift(n_nz:nframe,4)=(Drift(n_nz:nframe,4)-Drift(n_nz,4))*CCDpixelx_newton/CCDpixely;
%     Drift(n_nz:nframe,5)=Drift(n_nz:nframe,5)-Drift(n_nz,5);

    figure_drift = figure(1);
    set(figure_drift,  'Visible','Off');
    subplot(411), plot(1:nframe, Drift(:,3),'-'),title('dx0 (pixel)');
    subplot(412), plot(1:nframe, Drift(:,4),'-'),title('dy0 (pixel)');
    subplot(413), plot(1:nframe, Drift(:,5),'-'),title('dz0 (nm)');
    subplot(414), plot(1:nframe, Drift(:,1),'.'),title('dz0 (nm)');
    filePath_drift = [directoryName,'\Drift.fig'];
    saveas(figure_drift ,filePath_drift);
%
%
%
%Drift(:,1)=1; %??????????????????????????????????????????????????????
[n0,n]=size(find(Drift(:,1)==0));
[n1,n]=size(find(Drift(:,1)==1));
n2=n0+n1;
n=n1/n2*100;
fprintf('%g/%g (%2.0f percent) of images used for analysis. \n',n1, n2, n);
%
filePath_position = [directoryName,'\Drift.dat'];
fid_position = fopen(filePath_position,'w');
save(filePath_position,'Drift','-ASCII')
end