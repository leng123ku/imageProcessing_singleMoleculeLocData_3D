function q = curveFit2(p0)
global data_fit_input h

% fitWx = p0(1,1)+p0(1,2)*fitZ(:,1)+p0(1,3)*fitZ(:,2)+p0(1,4)*fitZ(:,3)+p0(1,5)*fitZ(:,4)+p0(1,6)*fitZ(:,5)+p0(1,7)*fitZ(:,6);
% fitWy = p0(2,1)+p0(2,2)*fitZ(:,1)+p0(2,3)*fitZ(:,2)+p0(2,4)*fitZ(:,3)+p0(2,5)*fitZ(:,4)+p0(2,6)*fitZ(:,5)+p0(2,7)*fitZ(:,6);
% 
% fitWx = p0(1,1)+p0(1,2)*fitZ(:,1)+p0(1,3)*fitZ(:,2)+p0(1,4)*fitZ(:,3)+p0(1,5)*fitZ(:,4);
% fitWy = p0(2,1)+p0(2,2)*fitZ(:,1)+p0(2,3)*fitZ(:,2)+p0(2,4)*fitZ(:,3)+p0(2,5)*fitZ(:,4);

% fitWx = p0(1,1)+p0(1,2)*fitZ(:,1)+p0(1,3)*fitZ(:,2)+p0(1,4)*fitZ(:,3);
% fitWy = p0(2,1)+p0(2,2)*fitZ(:,1)+p0(2,3)*fitZ(:,2)+p0(2,4)*fitZ(:,3);

fitZ = data_fit_input(:,1); 
fitWx = polyval(p0(1,:), data_fit_input(:,1));
fitWy = polyval(p0(2,:), data_fit_input(:,1));
[n,~]=size(data_fit_input);
q1 = 0;
q2 = 0;
% for i = 1:n
%     for j = 1:n
%         if fitZ(i,1) == data_fit_input(j,1)
%             q1=q1+sqrt((fitWx(i)-data_fit_input(j,2))^2+...
%                 (fitWy(i)-data_fit_input(j,3))^2);
%         end
%     end
% end

dist = sqrt((fitWx-data_fit_input(:,2)).^2 + (fitWy-data_fit_input(:,3)).^2);
q1 = sum(dist);


for i=1:n
    dist = (data_fit_input(i,2)-fitWx).^2+(data_fit_input(i,3)-fitWy).^2;
    minDis = min(dist);  
    q2 = q2 + minDis;
end
% k = q1/q2;
% fprintf('q1/q2 = %g\n', k);
q2 = 0*q2;

q = q1+q2;   % minimzing q
%q = q1/q2;

% if strcmpi(get(h.curveFit, 'Visible'), 'on')
%     figure(h.curveFit);
%     subplot(211); %Wx vs. Wy
%     plot(data_fit_input(:,2), data_fit_input(:,3), 'b.', ...
%         fitWx, fitWy, '.r');
%     legend('Data', 'Fit');
%     xlabel('Wx'); ylabel('Wy');
%     axis([1 3 1 3]); axis square
%     
%     subplot(212); %Wx vs. Wy
%     plot(data_fit_input(:,1), data_fit_input(:,2), 'b.', ...
%         data_fit_input(:,1), data_fit_input(:,3), 'g.', ...
%         fitZ(:,1), fitWx, '.r', ...
%         fitZ(:,1), fitWy, '.k');
%     %legend('Wx_Z_Data', 'Wy_Z_Data', 'Wx_Z_Fit', 'Wx_Z_Fit');
%     xlabel('Z'); ylabel('Wx or Wy');
%     axis([0 2000 1 5]); axis square
% end