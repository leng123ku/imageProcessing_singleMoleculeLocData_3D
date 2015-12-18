function [xc,yc,zc,dz,distance]  = cal_xyz14(x1,wx,dwx,y1,wy,dwy,ZWxy, p)
%% cal_xyz14 calculates Z of blinks (zc in nm), error in Z localization (dz in nm)
% and their corresponding distance (distance) using minimum distance method. 

Wx = ZWxy(:,1);
Wy = ZWxy(:,2);
Z = ZWxy(:,3);

% Finding minium distance to the calibration curve
[zc,distance] = distance01(wx, wy, Wx, Wy, Z);
% calculate dz
wx2 = wx;
wy2 = wy+dwy;
[zc1,~] = distance01(wx2, wy2, Wx, Wy, Z);

wx2 = wx;
wy2 = wy-dwy;
[zc2,~] = distance01(wx2, wy2, Wx, Wy, Z);

wx2 = wx+dwx;
wy2 = wy;
[zc3,~] = distance01(wx2, wy2, Wx, Wy, Z);

wx2 = wx-dwx;
wy2 = wy;
[zc4,~] = distance01(wx2, wy2, Wx, Wy, Z);

dz = (max(zc1,max(zc2,max(zc3,zc4)))-min(zc1,min(zc2,min(zc3,zc4))))./2;

%% Shift in x and y
%zc; 
% x_shift = polyval(px, zc);
% y_shift = polyval(py, zc);
% 
xc = x1;
yc = y1;


% x0=Coef(2,1)+Coef(2,2)+Coef(2,3)+Coef(2,4)+Coef(2,5)+Coef(2,6)+Coef(2,7);
% y0=Coef(3,1)+Coef(3,2)+Coef(3,3)+Coef(3,4)+Coef(3,5)+Coef(3,6)+Coef(3,7);
% Rxy=wx./wy;
% 
% x2=Coef(2,1)+Rxy.*(Coef(2,2)+ Rxy.*(Coef(2,3) + Rxy.*(Coef(2,4) + Rxy.*(Coef(2,5) + Rxy.*(Coef(2,6) + Rxy.*(Coef(2,7)))))));
% y2=Coef(3,1)+Rxy.*(Coef(3,2)+ Rxy.*(Coef(3,3) + Rxy.*(Coef(3,4) + Rxy.*(Coef(3,5) + Rxy.*(Coef(3,6) + Rxy.*(Coef(3,7)))))));
% 
% xc=x1-(x0-x2);
% yc=y1-(y0-y2);
% xc=x1;
% yc=y1;

