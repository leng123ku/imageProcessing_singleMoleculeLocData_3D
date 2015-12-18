function [position_fit]  = fit_blinks(directoryName)
%% fit_erf06 fits an error funtion to blinks

% - OUTPUT: 'position_fit'
%  position_fit(ndata, 18): fitted blinks info, used to contruct final image
%   description of columns
%   - 1:6 => bg, amp, x, wx, y, wy
%   - 7:12 => d_bd, d_amp, dx, dwx, dy, dwy
%   - 13 => R(goodness of fit)
%   - 14:15 => pixel(:,1:2)

% Written by Reza on June 17 2014
% Modified to obtain an efficient use of memory

%% Initialization
filePath_pixel = [directoryName,'\pixel.mat'];
px = matfile(filePath_pixel);
[n1, ~] = size(px.pixel);
[sizex, sizey, n2] = size(px.blinks);

if n1 ~= n2
    fprintf('Mismatch in size: pixel vs. blinks');
else
    nOfBlinks = n1;
    clear n1 n2
end

if sizex ~= sizey
    fprintf('Mismatch in cut size: not a square => check variable blinks');
else
    cutSize = sizex;
    clear sizex sizey
end

t_estimate = ceil(nOfBlinks/700/60);
t_now = clock;
t_now(4) = t_now(4) + floor( (t_now(5) + t_estimate) / 60);
t_now(5) = rem((t_now(5) + t_estimate) , 60);
cprintf('*red','Expecting to finish in %g min  =>  %g:%g. \n', t_estimate, t_now(4:5));

[Y, X] = meshgrid(1:cutSize);
x(:,1) = X(:);
x(:,2) = Y(:);

%% Fitting
noiselevel = 2;
%opts = optimset('TolFun',1e-3,'TolX',1e-3,'MaxIter', 30);
opts = optimset('TolFun',1e-6,'TolX',1e-6,'MaxIter', 100);
% fitFunc2 = @(c, surface) ...
%     c(1)+ c(2) *((erf((surface(:,1)-c(3)+0.5)/c(4))- ...
%     erf((surface(:,1)-c(3)-0.5)/c(4))) ...
%     .*(erf((surface(:,2)-c(5)+0.5)/c(6))-erf((surface(:,2)-c(5)-0.5)/c(6))));

fitFunc2 = @(c, surface) ...
    c(1) + ( c(2)/(c(4) * c(6) ) ) *exp(-( (surface(:,1) - c(3)).^2/ (2*c(4)^2) + ...
    ( surface(:,2) - c(5)).^2/ (2*c(6)^2) ));

% fitFunc2 = @(c, surface) ...
%     c(1)+c(2)*exp(-( (surface(:,1)-c(3))/c(4) ).^2 - ( (surface(:,2)-c(5))/c(6) ).^2);

%position = zeros(nOfBlinks,15);
position_fit = [];
filePath_position_fit = [directoryName,'\position_fit.mat'];
save(filePath_position_fit , 'position_fit');
n_FittedBlinks = 0;
warning off;

counter = 0;
blinkPerDivision = 5e4; % each dividion includes blinkPerDivision blinks
nOfDivision = ceil(nOfBlinks/blinkPerDivision);

t_totalFitting = tic;
fprintf('Currenly working on \n');
while counter < nOfDivision
    from = counter * blinkPerDivision + 1;
    if counter < nOfDivision - 1
        to = (counter+1) * blinkPerDivision;
    else
        to = nOfBlinks;
    end
    blinks = px.blinks(:,:, from:to);
    pixel = px.pixel(from:to,:);
    [n1, ~] = size(pixel);
    position_dev = zeros(n1,15);
    
    if counter == nOfDivision -1
            fprintf('%d to %d => LAST PORTION \n\n', from, to);
    else
        fprintf('%d to %d \n', from, to);
    end
    
%     if rem(counter,10) == 0
%         fprintf('%d to %d \n', from, to+9*blinkPerDivision);
%     else if counter == nOfDivision -1
%             fprintf('%d to %d => LAST PORTION \n', from, to);
%         end
%     end

    parfor k = 1:n1
        temp_position = zeros(1,15);
        temp_pixel = pixel(k,:);
        f = blinks(:,:,k)+1;
        [x0,y0,a,dx,dy] = gauss2dellipse(f,X,Y,noiselevel); % Anthony, Langmuir 2009, 25(14), 8152–8160
        ff = f(:);
        if (dx > 0.8) && (dy > 0.8)
            if (dx < 6) && (dy < 6)
                if ~((dx > 5) && (dy > 5))
                    a = a*dx*dy;
                    coeff = [3.5, a, x0, dx, y0, dy];
                    [fitCoeff, Res,~,CovB,~,~] = nlinfit( x, ff, fitFunc2, coeff, opts);
                    
                    conf2 = nlparci(fitCoeff,Res,'covar',CovB,'alpha',0.3173);
                    ci = conf2';
                    cfit = fitCoeff';
                    dc=(ci(2,:)-ci(1,:))/2;
                    totalRes = sum(( ff - mean(ff) ).^2);
                    R = 1 - sum(Res.^2)/totalRes;
                    
                    temp_position(1,1:6) = cfit;
                    temp_position(1,3) = cfit(3) - ceil(cutSize/2) + temp_pixel(1,3); % real x in actual image
                    temp_position(1,5) = cfit(5) - ceil(cutSize/2) + temp_pixel(1,4); % real y in actual image
                    temp_position(1,7:12) = dc;
                    temp_position(1,13) = R;
                    temp_position(1,14:15) = temp_pixel(1,1:2);
                    position_dev(k,:) = temp_position; % position_dev contains fitting info for the blinks in specific devision
                    n_FittedBlinks = n_FittedBlinks + 1;
                    
                end
            end
        end
    end
    position_dev = position_dev(any(position_dev,2),:); %cleaning rows with all zeros
    load(filePath_position_fit);
    position_fit = [position_fit; position_dev];
    save(filePath_position_fit , 'position_fit');
    position_fit = [];
    counter = counter + 1;
end

rate=n_FittedBlinks/toc(t_totalFitting);
timeInMin=toc(t_totalFitting)/60;
cprintf('*red', '%g/%g fitting finished in %.1f min at %g blinks/sec. \n\n',n_FittedBlinks, nOfBlinks, timeInMin, rate);

%% Writing to file, p0sition.dat
% position_fit = position_fit(any(position_fit,2),:); %cleaning rows with all zeros
% filePath_position = [directoryName,'\position_fit.dat'];
% save(filePath_position,'position_fit','-ASCII')

