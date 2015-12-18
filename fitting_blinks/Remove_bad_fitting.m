function [position] = Remove_bad_fitting(directoryName, nfile, nimage);
%%

%  position(ndata, 18): fitted blinks info, used to contruct final image
%   desc of columns
%   - 1:6 => bg, amp, x, wx, y, wy
%   - 7:12 => d_bd, d_amp, dx, dwx, dy, dwy
%   - 13 => R(goodness of fit)
%   - 14:15 => pixel(:,1:2)
%   - 16:18 => drift in x, y and z, respectively


%fprintf('Remove bad fittings (NaN) and add drift corrections...');
filePath= fullfile(directoryName,'\Drift.dat');
Drift = load(filePath);
% filePath= fullfile(directoryName,'\position_fit.dat');
% position_fit = load(filePath);
filePath= fullfile(directoryName,'\position_fit_falseSequentialBlinksRemoved.mat');
load(filePath);
position_fit = position_fit_falseSequentialBlinksRemoved;
[ndata, np]=size(position_fit);
clear position_fit_falseSequentialBlinksRemoved
position=zeros(ndata, np+3);
point=0;
for k=1:ndata
    if (position_fit(k,14) ~= 0)
        [n1,n2]=size(find(position_fit(k,:)>=0 & position_fit(k,:) < 1e8));
        if n2 == 15
            point=point+1;
            position(point,1:np)=position_fit(k,:);
            [n_newton, n2]=size(Drift);
            Nf=position_fit(k,14); %file No.
            nD=sum(nimage(1:Nf-1))+position_fit(k,15);  % image No. in Newton
%            position(point,np+1:np+3)=Drift(nD,3:5);
        end 
    end
end
if point == 0
    fprintf('No blinks!')
end
position = position(any(position,2),:);
%position_out=position(1:point,:);
% filePath_position = [directoryName,'\position.dat'];
% fid_position = fopen(filePath_position,'w');
% save(filePath_position,'position_out','-ASCII');
filePath_position = [directoryName,'\position.mat'];
save(filePath_position,'position');
removed=ndata-point;
cprintf('*red', '%g blinks removed => %g points left. \n\n', removed, point);
fclose('all');