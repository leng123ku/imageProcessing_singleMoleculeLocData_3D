% calculate z position using Rxy and shift x and y according to the z
% position
function [xc,yc,zc]  = cal_xyz03(x1,wx,y1,wy,Coef)
%
%
% x0=Coef(2,1)+Coef(2,2)+Coef(2,3)+Coef(2,4)+Coef(2,5)+Coef(2,6)+Coef(2,7);
% y0=Coef(3,1)+Coef(3,2)+Coef(3,3)+Coef(3,4)+Coef(3,5)+Coef(3,6)+Coef(3,7);
Rxy=wx./wy;
zc=Coef(1,1)+Rxy.*(Coef(1,2)+ Rxy.*(Coef(1,3) + Rxy.*(Coef(1,4) + Rxy.*(Coef(1,5) + Rxy.*(Coef(1,6) + Rxy.*(Coef(1,7)))))));
% x2=Coef(2,1)+Rxy.*(Coef(2,2)+ Rxy.*(Coef(2,3) + Rxy.*(Coef(2,4) + Rxy.*(Coef(2,5) + Rxy.*(Coef(2,6) + Rxy.*(Coef(2,7)))))));
% y2=Coef(3,1)+Rxy.*(Coef(3,2)+ Rxy.*(Coef(3,3) + Rxy.*(Coef(3,4) + Rxy.*(Coef(3,5) + Rxy.*(Coef(3,6) + Rxy.*(Coef(3,7)))))));
%
% xc=x1-(x0-x2);
% yc=y1-(y0-y2);
xc=x1;
yc=y1;