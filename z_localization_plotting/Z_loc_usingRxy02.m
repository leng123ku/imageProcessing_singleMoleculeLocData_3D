function [zc, dz]  = Z_loc_usingRxy02(Wx,dWx,Wy,dWy, Rxy_Z_table)
%% cal_xyz06: calculates zc and dz using minimum distance method. 
% - xy shift is included. 
% - should be used in conjuction with plot_blinks04 and storm_new06
% beacause of Coef. 

Rxy_min = min(Rxy_Z_table(:,1));
Rxy_max = max(Rxy_Z_table(:,1));
% calculate zc
Rc = Wx./Wy;
[n, ~] = size(Rc);
zc = zeros(n,1);

for i = 1:n
    if ((Rc(i) > Rxy_min) && (Rc(i) < Rxy_max))
        
        R_dif = abs(Rc(i) - Rxy_Z_table(:,1));
        [~, m] = min(R_dif);
        zc(i,1) = Rxy_Z_table(m,2);
        
    else
        zc(i,1) = NaN;
    end
end

% calculate dz
Wy_up= Wy + dWy;
Rc = Wx./Wy_up;
zc_Wy_up = zeros(n,1);

for i = 1:n
    if ((Rc(i) > Rxy_min) && (Rc(i) < Rxy_max))
        
        R_dif = abs(Rc(i) - Rxy_Z_table(:,1));
        [~, m] = min(R_dif);
        zc_Wy_up(i,1) = Rxy_Z_table(m,2);
        
    else
        zc_Wy_up(i,1) = NaN;
    end
end


Wy_dn= Wy - dWy;
Rc = Wx./Wy_dn;
zc_Wy_dn = zeros(n,1);

for i = 1:n
    if ((Rc(i) > Rxy_min) && (Rc(i) < Rxy_max))
        
        R_dif = abs(Rc(i) - Rxy_Z_table(:,1));
        [~, m] = min(R_dif);
        zc_Wy_dn(i,1) = Rxy_Z_table(m,2);
        
    else
        zc_Wy_dn(i,1) = NaN;
    end
end


Wx_up= Wx + dWx;
Rc = Wx_up./Wy;
zc_Wx_up = zeros(n,1);

for i = 1:n
    if ((Rc(i) > Rxy_min) && (Rc(i) < Rxy_max))
        
        R_dif = abs(Rc(i) - Rxy_Z_table(:,1));
        [~, m] = min(R_dif);
        zc_Wx_up(i,1) = Rxy_Z_table(m,2);
        
    else
        zc_Wx_up(i,1) = NaN;
    end
end

Wx_dn= Wx - dWx;
Rc = Wx_dn./Wy;
zc_Wx_dn = zeros(n,1);

for i = 1:n
    if ((Rc(i) > Rxy_min) && (Rc(i) < Rxy_max))
        
        R_dif = abs(Rc(i) - Rxy_Z_table(:,1));
        [~, m] = min(R_dif);
        zc_Wx_dn(i,1) = Rxy_Z_table(m,2);
        
    else
        zc_Wx_dn(i,1) = NaN;
    end
end


dz=(max(zc_Wy_up,max(zc_Wy_dn,max(zc_Wx_up,zc_Wx_dn)))-min(zc_Wy_up,min(zc_Wy_dn,min(zc_Wx_up,zc_Wx_dn))))./6;
% xc = x1;
% yc = y1;


% x_shift = polyval(px, zc);
% y_shift = polyval(py, zc);
% 
% xc = x1 - x_shift;
% yc = y1 - y_shift;


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

