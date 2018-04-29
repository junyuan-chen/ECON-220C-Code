tic
clear;
close all;
N = 100;
T = 6;
nreplic = 1000;
rho0 = 0.7;


NT=N*T;
 

% keep the bias of each estimator
bias_fe=zeros(10,1);
bias_ah1=zeros(10,1);
bias_ah2=zeros(10,1);  

% keep the se of each estimator
se_fe = zeros(10,1);
se_ah1 = zeros(10,1);
se_ah2 = zeros(10,1);

% keep the rmse of each estimator
rmse_fe = zeros(10,1);
rmse_ah1 = zeros(10,1);
rmse_ah2 = zeros(10,1);


rho = (0.1:0.1:1)'; %the values of rho considered

for irho = 1:1:10;  % loop over each value of rho's

  a = [1 -rho(irho)];
  b = 1;  % configure a and b to speed up the data generating

  rho_fe = zeros(nreplic,1);  % keep the fixed effects estimator
  rho_ah1 = zeros(nreplic,1);
  rho_ah2 = zeros(nreplic,1);
  
  t_fe = zeros(nreplic,1); % keep the t statistic based on the fe estimator
  t_ah1 = zeros(nreplic,1);
  t_ah2 = zeros(nreplic,1);
  
  randn('state',pi);
% initialize the random number generator

   for replic = 1:1:nreplic; 
    
		alpha = randn(1,N);  
		e = randn(1,N);
		y0 = rho0*alpha+e; 
		
	
		temp = [y0; randn(T,N)+alpha(ones(T,1),:)];
		data  =  filter(b,a, temp);   
		
		x = data(1:end-1,:);
		y = data(2:end,:);
          
		% data in wide format
          
		clear temp data;
         

			
		%%%%%%%%%%%%% FE estimator %%%%%%%%%%%%%%%%%
		x_fe = x-repmat(mean(x),T,1);
		y_fe = y-repmat(mean(y),T,1);
		rho_fe(replic) = x_fe(:)\y_fe(:);
        
        
        err_fe = y_fe(:)-x_fe(:)*rho_fe(replic); % regression error
            
        
        avar_fe = 0;
        for ii=1:1:N;
         
        index = ((ii-1)*T+1: ii*T)';
         avar_fe = avar_fe + x_fe(index)'*err_fe(index)*err_fe(index)'*x_fe(index);
        end;
        
        avar_fe = inv(x_fe(:)'*x_fe(:))*avar_fe*inv(x_fe(:)'*x_fe(:));
        
         t_fe(replic) = (rho_fe(replic)-rho(irho))/sqrt(avar_fe);
		

		
		%%%%%%%%%%%%% AH1 estimator %%%%%%%%%%%%%%%%%

        
        x_fd = x(2:end,:)-x(1:end-1,:);
		y_fd = y(2:end,:)-y(1:end-1,:);
        
        
		z = reshape(x(1:end-1,:),[],1);


        rho_ah1(replic) = (z'*y_fd(:))/(z'*x_fd(:));
        
        %     'Need to compute  t_ah1(replic)' 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        		
		%%%%%%%%%%%%% AH2 estimator %%%%%%%%%%%%%%%%%
        
        x_fd = x(3:end,:)-x(2:end-1,:);
		y_fd = y(3:end,:)-y(2:end-1,:);
        
        z = reshape(x(2:end-1,:)-x(1:end-2,:),[],1);
        
        
		rho_ah2(replic) = (z'*y_fd(:))/(z'*x_fd(:));
        
             % 'Need to compute  t_ah2(replic)' 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
    end;

bias_fe(irho)  =  mean(rho_fe)-rho(irho);
bias_ah1(irho)  =  mean(rho_ah1)-rho(irho);
bias_ah2(irho)  =  mean(rho_ah2)-rho(irho);


se_fe(irho)  =  std(rho_fe);
se_ah1(irho)  =  std(rho_ah1);
se_ah2(irho)  =  std(rho_ah2);


rmse_fe(irho)  =  sqrt(se_fe(irho)^2+bias_fe(irho)^2);
rmse_ah1(irho)  =  sqrt(se_ah1(irho)^2+bias_ah1(irho)^2);
rmse_ah2(irho)  =  sqrt(se_ah2(irho)^2+bias_ah2(irho)^2);

	
	if rho(irho) == 0.7;
		figure(irho);
   
		[freq, bins]= hist(t_fe);
        bar(bins, freq/sum(freq));
        hold on;
        x=(-2.5:0.1:2.5)';
        plot(x, normpdf(x,0,1));
		title({'Histogram of the t_{FE} statistic'; ['\rho = ' num2str(0.7) ', N = ' num2str(N) ', T = ' num2str(T)]});
   
        orient landscape;
        print('-depsc2', ['histogram' '_N_' num2str(N) '_T_' num2str(T) '.eps'])
      
        % 'Need to graph the histogram for other estimators'
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %                                                       %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end;
     

end;

figure(100)

H1 = plot(rho, bias_fe, 'k+:',...
     rho, bias_ah1,  'bo-',...
     rho, bias_ah2,  'rs-.');
H2= legend('fe','ah1','ah2');
title(['Bias of Different Estimators,' ' N = ' num2str(N) ', T = ' num2str(T)]);
xlabel ('\rho'); 
ylabel ('bias');
set(H1, 'linewidth', 2,'markersize', 4);
set(H2, 'fontsize', 12)


orient landscape;
print('-depsc2', ['bias_se_rmse' '_N_' num2str(N) '_T_' num2str(T) '.eps'])

        %'Need to draw the graphs for other estimators'
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %                                                       %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       
toc
