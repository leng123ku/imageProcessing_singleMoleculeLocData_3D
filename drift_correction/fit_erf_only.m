function [P]  = fit_erf_only(image)
%
%
[nsize,nsize2,ndata]=size(image);
%
%nshow=50;  %<===============================
%nd=100;
P=zeros(ndata,14);
%image_temp=zeros(nsize,nsize2);
%
for n=1:nsize
    for m=1:nsize2
        X(n,m)=n;
        Y(n,m)=m;
    end
end
%[X,Y]=meshgrid(1:nsize,1:nsize);%x-y coordinates
x(:,1)=X(:); % x= first column
x(:,2)=Y(:); % y= second column
%
fun = fittype( 'a+b*((erf((x-c+0.5)/d)-erf((x-c-0.5)/d)).*(erf((y-e+0.5)/f)-erf((y-e-0.5)/f)))',...
      'coefficients', {'a', 'b', 'c', 'd', 'e', 'f'},'independent', {'x', 'y'},'dependent', 'z' );
% Set initial fitting parameters
image1_max=max(max(image(:,:,1)));
% [x0,y0]=find(image(:,:,1) == image1_max);
% dx=1.7;            % dx
% dy=1.7;            % dy
% a0=1;
tic
% fitting +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
parfor k=1:ndata  
    image_temp=image(:,:,k)+1;
    f=image_temp(:);          % your data f(x,y) (in column vector)
    temp_position=zeros(1,14);
    %
    noiselevel=min(f(:))*2;
    [x0,y0,a,dx,dy]=gauss2dellipse(image_temp,X,Y,noiselevel); % Anthony, Langmuir 2009, 25(14), 8152–8160
    a=a*3.6724;
    a0=noiselevel/3;
    FlagNaN=sum(isnan([a0, a, x0, dx, y0, dy])); % Check the estimate
    if FlagNaN ~= 0         % if NaN appear, then do not fit
        %fprintf('image No. %g contains NaN \n',k)
        temp_position(1,1:14)=NaN;
        P(k,:)=temp_position;
    else 
        %  
        [cc,gof,output] = fit(x,f,fun,...
                   'Lower',[0,1,1,0.5,1,0.5],...
                   'Upper',[Inf,Inf,nsize,5,nsize,5],...
                   'Startpoint',[a0, a, x0, dx, y0, dy],...
                   'Display','off',...
                   'MaxIter',1e2,... %100
                   'TolFun',1e-10,... %1e-2
                   'TolX',1e-4);     %1e-2
        cfit=coeffvalues(cc); 
        ci = confint(cc,0.683);
        dc=(ci(2,:)-ci(1,:))/2;
        dc(4)=dc(4)/sqrt(2);
        dc(6)=dc(6)/sqrt(2);
    %
        % record positions  ---------------------------------
        % Check exitflag
        %iteration=output.iterations;
        exitflag=output.exitflag; %Positive flags indicate convergence, within tolerances.
        exit_iteration=output.iterations;
        R=gof.rsquare;  %Coefficient of determination
        if R < 0.9 | R > 1.0 | cfit(2) < 10
            %fprintf('image No. %g give bad fit R or low intensity\n',k)
            temp_position(1,1:14)=NaN;
            P(k,:)=temp_position;
        elseif sum(isnan(cfit)) ~= 0 | sum(isnan(dc)) ~= 0
            %fprintf('image No. %g give bad fitting parameters \n',k)
            temp_position(1,1:14)=NaN;
            P(k,:)=temp_position;
        else
            temp_position(1,1:6)=cfit;
            temp_position(1,4)=cfit(4)/sqrt(2);
            temp_position(1,6)=cfit(6)/sqrt(2);
            temp_position(1,7:12)=dc;
            temp_position(1,13)=R;
            temp_position(1,14)=k;
            P(k,:)=temp_position;
        end
           %
    %      if k <= 1
%     fprintf('k: %g \n',k)
%             fprintf('Initial parameter: [%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f] \n',a0, a, x0, dx, y0, dy)
%             fprintf('Fianl   parameter: [%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f] \n',cfit)
%             fprintf('Fianl        dc  : [%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f] \n',dc)
%              fprintf('exit_iteration = %g; exitflag= %g, ; R= %g \n',exit_iteration, exitflag, R)
%             % Plot
%             A=zeros(nsize,nsize2);
%             I=cfit(2);
%             bg=cfit(1);            
%             Xpos=cfit(3);
%             Ypos=cfit(5);
%             sigmax=cfit(4);%2+(n-1)*0.4;
%             sigmay=cfit(6);%6-wx;
%              A=I*((erf((X-Xpos+.5)./(sqrt(2)*sigmax))-erf((X-Xpos-.5)./(sqrt(2)*sigmax))).*...
%         (erf((Y-Ypos+.5)./(sqrt(2)*sigmay))-erf((Y-Ypos-.5)./(sqrt(2)*sigmay))))+bg;
%             figure(1) 
%             subplot(211), imshow(A,[min(A(:)) max(A(:))]), title('fit');
%             subplot(212), imshow(image(:,:,k),[min(A(:)) max(A(:))]), title('data');
    %      end
    end
end %============================================
time=toc;
rate=ndata/time;
fprintf('Finished at %g fit/sec\n',rate);