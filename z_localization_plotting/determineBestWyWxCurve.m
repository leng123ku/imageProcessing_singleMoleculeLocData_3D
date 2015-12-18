function [Wx_Wy_Z] = determineBestWyWxCurve(selBlinks, Rxy_Z_table, directoryName, ZWxy)
global data_fit_input

[n1, ~] = size(selBlinks);
beadRemoval = false(n1,1);
Rc = selBlinks(:,1)./selBlinks(:,2);

% Rxy_min = min(Rxy_Z_table(:,1));
% Rxy_max = max(Rxy_Z_table(:,1));
% 
% 
% for i = 1:n1
%     if ((Rc(i) > Rxy_min) && (Rc(i) < Rxy_max))
%         
%         R_dif = abs( Rc(i) - Rxy_Z_table(:,1) );
%         [~, m] = min(R_dif);
%         selBlinks(i,3) = Rxy_Z_table(m,2);
%     else
%         selBlinks(i,3) = NaN;
%         beadRemoval(i) = true;
%     end
% end


ulimWx = max(ZWxy(:,1));
ulimWy = max(ZWxy(:,2));

for i = 1:n1
    if ((selBlinks(i,1) < ulimWx) && (selBlinks(i,2) < ulimWy))
        
        R_dif = abs( Rc(i) - Rxy_Z_table(:,1) );
        [~, m] = min(R_dif);
        selBlinks(i,3) = Rxy_Z_table(m,2);
    else
        selBlinks(i,3) = NaN;
        beadRemoval(i) = true;
    end
end

toRemove = find(beadRemoval == true);
selBlinks(toRemove,:) = [];
Rc(toRemove,:) = [];

[n2, ~] = size(selBlinks);
fprintf('%d / %d blinks used for finding an optimum calibration curve\n', n2, n1);


h1 = figure; plot(selBlinks(:,3), Rc, '.b');
ylabel('Z'); xlabel('Rxy');
fid = [directoryName '\Rxy_Z_blinks.fig'];
saveas(h1, fid);

[n, ~] = size(selBlinks);
data_fit_input = zeros(n,3);
data_fit_input(:,1) = selBlinks(:,3);
data_fit_input(:,2) = selBlinks(:,1);
data_fit_input(:,3) = selBlinks(:,2);

minZ = min(data_fit_input(:,1));
maxZ = max(data_fit_input(:,1));

% fittedData: Z, Wx and Wy of beads, fitting results
% - (:,1): Z, (:,2): Wx, (:,3): Wy
Z = (minZ:0.5:maxZ)';  % scliced Z to get wx_fit and wy_fit
%Z = (90:1000)';
[n, ~] = size(Z);
data_fit_output = zeros(n,3);
data_fit_output(:,1) = Z;


% p0: initial coeff for fitting a 3rd-degree polynomial
polyDeg = 3;
p1 = polyfit(data_fit_input(:,1), data_fit_input(:,2), polyDeg);
%p1 = fliplr(p1); % flip the coefficients/mir

p2 = polyfit(data_fit_input(:,1), data_fit_input(:,3), polyDeg);
%p2 = fliplr(p2);

p0=0;
p0(1,1:polyDeg+1) = p1;
p0(2,1:polyDeg+1) = p2;

% data_fit_output(:,2) = polyval(p0(1,:), Z);
% data_fit_output(:,3) = polyval(p0(2,:), Z);

options = optimset('MaxFunEvals',2e6,'MaxIter',2e3,'TolFun',1e-2,'TolX',1e-2,'display','final');
%options = optimset('MaxIter',1e5,'TolFun',1e-15,'TolX',1e-15,'display','final');

h.curveFit = figure('Visible','off');
for i=1 % i: number of times to repeat the curve fitting
    [pf,q,exitflag,output] = fminsearch(@(x) curveFit4(x),p0, options);
    p0=pf;
    pf = p0;
    %fprintf('q= %g; q1= %g; q2= %g \n',q, q1, q2);
    fprintf('     p0(1,:)=[%E %E %E %E %E %E %E]; \n',pf(1,:));
    fprintf('     p0(2,:)=[%E %E %E %E %E %E %E];\n',pf(2,:));
    %fprintf(' exitflag = %g -- iterations = %g\n', exitflag, output.iterations);
    data_fit_output(:,2) = polyval(pf(1,:), Z);
    data_fit_output(:,3) = polyval(pf(2,:), Z);

end

% Writing calibration file, ZWxy => Wx, Wy, Z.
fprintf('Writing calibration file\n');
[n, ~] = size(data_fit_output);
ZWxy1 = zeros(n,4);
k = 1;
for i=1:n
    if (data_fit_output(i,2) <=3) && (data_fit_output(i,3) <=3) % write only when Wx,Wy<=3
        
        ZWxy1(k,1) = data_fit_output(i,2);
        ZWxy1(k,2) = data_fit_output(i,3);
        ZWxy1(k,3) = data_fit_output(i,1);
        ZWxy1(k,4) = ZWxy1(k,1)./ZWxy1(k,2);
        k = k+1;
    end
end
ZWxy1 = ZWxy1(any(ZWxy1,2),:); % remove row with all zeros
Wx_Wy_Z = ZWxy1;

h1 = figure;
plot(data_fit_input(:,2), data_fit_input(:,3), 'b.', ...
    data_fit_output(:,2), data_fit_output(:,3), '.r', Wx_Wy_Z(:,1), Wx_Wy_Z(:,2), 'g.');
h_legend = legend('Data', 'Fit', 'Selection-ZWxy');
set(h_legend,'FontSize',12);
xlabel('Wx'); ylabel('Wy');
axis([0 4 0 4]); axis square

fid = [directoryName '\Wy_Wx_blinks_fit.fig'];
saveas(h1, fid);

h1 = figure;
xedges = linspace(0,4,100); yedges = linspace(0,4,100);
histmat = hist2(data_fit_input(:,2), data_fit_input(:,3), xedges, yedges);
pcolor(xedges,yedges,histmat'); colorbar ; axis square tight; 
hold on
plot(data_fit_output(:,2), data_fit_output(:,3), '.r', Wx_Wy_Z(:,1), Wx_Wy_Z(:,2), 'g.');
xlabel('Wx'); ylabel('Wy');
axis([0 4 0 4]); axis square
fid = [directoryName '\Wy_Wx_blinks_fit_hist.fig'];
saveas(h1, fid);


h1 = figure;
plot(data_fit_input(:,1), data_fit_input(:,3), 'b.', ...
    data_fit_input(:,1), data_fit_input(:,2), 'g.', ...
    data_fit_output(:,1), data_fit_output(:,3), '.r', ...
    data_fit_output(:,1), data_fit_output(:,2), '.k');
h_legend = legend('Wy-Data', 'Wx-Data', 'Wy-Fit', 'Wx-Fit');
set(h_legend,'FontSize',12);
xlabel('Z'); ylabel('Wx or Wy');
axis([-1000 2000 0.5 5]); axis square

fid = [directoryName '\Wy_Z_Wx_Z_fitPolynomial.fig'];
saveas(h1, fid);














