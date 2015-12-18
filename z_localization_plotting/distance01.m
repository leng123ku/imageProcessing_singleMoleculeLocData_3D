function [zc,distance]  = distance01(wx, wy, Wx, Wy, Z)
%% distance01 calculates Z of blinks and their corresponding distance to calibration curve
% INPUT: 
% - wx, wy: width of blinks (vector)
% - Wx, Wy, Z: calibration curve info (vector)
% OUTPUT:
% - zc: Z of blink (vector)
% - distance : corresponding distance to calibration curve (vector)


[nOfBlinks,~] = size(wx);
distance = zeros(nOfBlinks,1);
zc = zeros(nOfBlinks,1);

parfor n = 1:nOfBlinks
    ZZ = Z;
    %D = sqrt((wx(n,1)-Wx).^2+(wy(n,1)-Wy).^2);
    D = sqrt((wx(n,1).^(0.5)-Wx.^(0.5)).^2+(wy(n,1).^(0.5)-Wy.^(0.5)).^2);
    minD = min(D);
    [row_index,~] = find(D == minD);
    [noOfMinDis, ~] = size(row_index);
    if noOfMinDis == 1
        zc(n,1) = ZZ(row_index,1);
        distance(n,1) = minD;
    else
        zc(n,1) = NaN;
        distance(n,1) = NaN;
    end
end
