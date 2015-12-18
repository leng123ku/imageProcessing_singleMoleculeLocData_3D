function cut  = cut_image(directoryName, ntx, nty, cut_edge)
%
preSTORM_path = fullfile(directoryName,'PreSTORM.tif'); %< input wide-filed Images +++++++
flag_pre = exist(preSTORM_path, 'file');
if flag_pre ~= 0
    I = (preSTORM_path);
    info = imfinfo(I);
    %num_images = numel(info);
    Pre_STORM = imread(I, 1, 'Info', info);
    STORM_1 = sum(Pre_STORM,1); % 1x512
    STORM_2 = sum(Pre_STORM,2); % 512x1
    dSTORM_1 = diff(STORM_1); % 1x512
    dSTORM_2 = diff(STORM_2); % 512x1
    %         DSTORM_TH = min([mean(abs(dSTORM_1(1,1:10))) mean(abs(dSTORM_1(1,nty-11:nty-1)))...
    %             mean(abs(dSTORM_2(1:10,1))) mean(abs(dSTORM_2(ntx-11:ntx-1,1))) ])*10;
    cut(1)=1;  % Top (pixel)  Y
    cut(2)=ntx-1;  % buttom (pixel) Y
    cut(3)=1;  % Left (pixel) X
    cut(4)=nty-1;  %Right (pixel) X
    STORM_TH1 = max(STORM_1)/12;
    STORM_TH2 = max(STORM_2)/12;
    while STORM_2(cut(1),1) < STORM_TH2
        cut(1)=cut(1)+1;
    end
    cut(1) = cut(1)-7;
    if cut(1)<1
        cut(1)=1;
    end
    while STORM_2(cut(2),1) < STORM_TH2
        cut(2) = cut(2)-1;
    end
    cut(2) = cut(2)+7;
    if cut(2) > ntx
        cut(2) = ntx;
    end
    while STORM_1(1,cut(3)) < STORM_TH1
        cut(3)=cut(3)+1;
    end
    cut(3) = cut(3)-7;
    if cut(3) < 1
        cut(3) = 1;
    end
    while STORM_1(1,cut(4)) < STORM_TH1
        cut(4)=cut(4)-1;
    end
    cut(4) = cut(4)+7;
    if cut(4) > nty
        cut(4)=nty;
    end
    %         Pre_STORM_cut = Pre_STORM(cut(1):cut(2),cut(3):cut(4));
    
    %         figure(100)
    %         subplot(321),imshow(Pre_STORM,[min(Pre_STORM(:)) max(Pre_STORM(:))]), title('Pre_STORM');
    %         subplot(322), plot(1:ntx,STORM_1,'.-'), title('STORM_1');
    %         subplot(323), plot(1:nty,STORM_2,'.-'), title('STORM_2');
    %         subplot(324), plot(1:ntx-1,dSTORM_1,'.-'), title('dSTORM_1');
    %         subplot(325), plot(1:nty-1,dSTORM_2,'.-'), title('dSTORM_2');
    %         subplot(326),imshow(Pre_STORM_cut,[min(Pre_STORM_cut(:)) max(Pre_STORM_cut(:))]), title('Pre_STORM_cut');
elseif flag_pre == 0
    cutx = 0;%ntx-20*floor(ntx/20);
    cuty = 0;%nty-20*floor(nty/20);
    cut(1) = 1;  % Top (pixel)  Y
    cut(3) = 1;  % Left (pixel) X
    cut(2) = ntx-cutx;  % buttom (pixel) Y
    cut(4) = nty-cuty;  %Right (pixel) X
end
%
cut(1) = 1;%max(cut(1), cut_edge);
cut(2) = 512;%min(cut(2), 512-cut_edge+20);
cut(3) = 1;%max(cut(3), cut_edge);
cut(4) = 512;%min(cut(4), 512-cut_edge+20);
% 
% cut(1) = max(cut(1), cut_edge)
% cut(2) = min(cut(2), 512-cut_edge+20)
% cut(3) = max(cut(3), cut_edge)
% cut(4) = min(cut(4), 512-cut_edge+20)

%% writing the widefield image.
preSTORM_path = fullfile(directoryName,'PreSTORM.tif'); %< input wide-filed Images +++++++
flag_pre = exist(preSTORM_path, 'file');
if flag_pre == 2
    info = imfinfo(preSTORM_path);
    A = imread(preSTORM_path, 1, 'Info', info);
    A_cut = A(cut(1):cut(2),cut(3):cut(4));
    wideField_path = fullfile(directoryName,'WideField.tif');
    imwrite(A_cut,wideField_path,'Compression','none');
else
    error('No pre_STORM image was found in %s', preSTORM_path);
end




