clear;
clc;
N = 500;
T = 5;
nreplic = 1000;

beta_hat = zeros(nreplic,1);
sigma_beta1 =  zeros(nreplic,1);
sigma_beta2 =  zeros(nreplic,1);

for replic = 1:1:nreplic;
    
    x = randn(T,N);
    %sigma = repmat(sqrt(T)*mean(x),T,1);
    sigma = 2*abs(x);
    e = sigma.*randn(T,N);
    y = x + e;
    
    
    y = y - repmat(mean(y),T,1);
    x = x - repmat(mean(x),T,1);
    
    yy = y';  % T x N
    xx = x';  % T x N
    
    xx = xx(:); % vectorize x' 
    yy = yy(:); % vectorize y'

    
    beta_hat(replic) = xx\yy;  %fixed effects estimator
    e_hat = y - x*beta_hat(replic);  % residual from fixed effect regression
                                     % T x N
  
    
    sigma_beta1(replic) = sqrt((xx'*xx)^(-2)*sum( sum(x.*e_hat).^2));
    sigma_beta2(replic) = sqrt((xx'*xx)^(-2)*sum(sum ( (x.^2).*(e_hat.^2) )));
end


[std(beta_hat) mean(sigma_beta1) mean(sigma_beta2) std(sigma_beta1), std(sigma_beta2)...
    sqrt((mean(sigma_beta1)-std(beta_hat))^2 +std(sigma_beta1)^2), sqrt((mean(sigma_beta2)-std(beta_hat))^2 +std(sigma_beta2)^2)   ]


figure(10)
hist(beta_hat);  

figure(1)
hist(sigma_beta1);

figure(2)
hist(sigma_beta2);
    

