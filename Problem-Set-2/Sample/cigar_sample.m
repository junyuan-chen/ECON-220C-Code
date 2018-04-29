%% Load the data
clear;
clc;
cigar_new =load('cigar_new.csv');

state=cigar_new(:,1);	
yr=cigar_new(:,2);
lnc=cigar_new(:,3);
lnp=cigar_new(:,4);

lny=cigar_new(:,5);
lnpn=cigar_new(:,6);

cig=[lnp lny lnpn yr lnc];
[rows, cols]=size(cig); 

%% Mean and Demean 
T=30; % length of time series
N=rows/T; % number of cross-sectional units
a=ones(T-1,1); 

data = [];
data_dm = []; %variable for demeaned data
data_m = []; %variable for cross sectional mean

for i = 1:1:N;
    data_temp = [cig((i-1)*T+1:(i*T)-1,cols) cig((i-1)*T+2:i*T,:)];
    data = [data; data_temp];
    
    data_m = [data_m; a*mean(data_temp)];
    data_dm = [data_dm; data_temp-a*mean(data_temp)];
end;

T = T-1; 
[~, cols]=size(data);

%% pooled OLS estimate with no year dummies 

y = data(:,cols);
x = [data(:,1:4) ones(length(y),1)];

beta_ols = x\y;
u_ols = y-x*beta_ols;
[~,K] = size(x);
var_beta_ols =(u_ols'*u_ols)/(N*T-5)*((x'*x)\eye(K,K));
se_beta_ols = sqrt(diag(var_beta_ols));
t_beta_ols = beta_ols./se_beta_ols;



%%  pooled OLS with year dummy 
yr=data(:,cols-1);
unique_yr=unique(yr);
nyr=length(unique_yr); % number of years in the data;
dummy=[];
for i = 2:1:nyr;  %to avoid the dummy varible trap, i starts from 2
    dummy = [dummy (yr==unique_yr(i))];
end; 
                  % the loop here is not the optimal coding
                  % the size of 'dummy' grows with i

x=[data(:,1:4) dummy ones(length(y),1)];
beta_ols_d=x\y;
[~,K] = size(x);
u_ols_d=y-x*beta_ols_d;
var_beta_ols_d=(u_ols_d'*u_ols_d)/(N*T-4-nyr)*((x'*x)\eye(K,K));
se_beta_ols_d=sqrt(diag(var_beta_ols_d));
t_beta_ols_d=beta_ols_d./se_beta_ols_d;



%%  within estimate without year dummies 
y=data_dm(:,cols);
x=data_dm(:,1:4);

beta_within=x\y;
[~,K] = size(x);
u_within=y-x*beta_within;
var_beta_within=(u_within'*u_within)/(N*(T-1)-4)*((x'*x)\eye(K,K));
se_beta_within=sqrt(diag(var_beta_within));
t_beta_within=beta_within./se_beta_within;



%%  within estimate with year dummies 
y=data_dm(:,cols);
x=[x ones(length(y),1) dummy];

beta_within_d=x\y;
[~,K] = size(x);
u_within_d=y-x*beta_within_d;
var_beta_within_d=(u_within_d'*u_within_d)/(N*(T-1)-4-nyr)*((x'*x)\eye(K,K));
se_beta_within_d=sqrt(diag(var_beta_within_d));
t_beta_within_d=beta_within_d./se_beta_within_d;


%% Anderson and Hsiao Estimator 
[~, cols]=size(cig);

data_fd=[];
iv=[];
T=30;
for i=1:1:N;
    y_temp3=cig((i-1)*T+3:(i*T),cols); %3:T
    y_temp2=cig((i-1)*T+2:(i*T)-1,cols); %2:T-1
    y_temp1=cig((i-1)*T+1:(i*T)-2,cols); %1:T-2
    delta_y3=y_temp3-y_temp2;
    delta_y2=y_temp2-y_temp1;
    
    data_temp=[delta_y2 cig((i-1)*T+3:i*T,1:3)-cig((i-1)*T+2:(i*T)-1,1:3) delta_y3];
    
    data_fd=[data_fd; data_temp];
    iv=[iv; y_temp1];
  
end;

T=29;
dummy2=[]; 
dummy3=[];

for i=1:1:N;
  dummy2=[dummy2; dummy((i-1)*T+2:(i*T),:)-dummy((i-1)*T+1:(i*T-1),:)];
  dummy3=[dummy3; dummy((i-1)*T+2:(i*T),:)];
end;

[~, cols]=size(data_fd);
y=data_fd(:,cols);
x=[data_fd(:,1:4) dummy3]; %it does not matter whether dummy2 or dummy 3 is used. 
z=[iv data_fd(:,2:4) dummy3];

%Pz=z*inv(z'*z)*z';
%beta_ah=inv(x'*Pz*x)*(x'*Pz*y); %exactly identified
beta_ah=(z'*x)\(z'*y);



%% Arellano and Bond Estimator using all the instruments 
% first construct the instrument matrix
% using the lagged lnc up to t=1

Z=[];
Zi=[];

T=30; 
for i=1:1:N;
 yi=lnc((i-1)*T+1:i*T,:);
 Zi=[];
   for t=1:1:T-2
       Zi=blkdiag(Zi,yi(1:t)');
    end;
Z=[Z;Zi];
end;

% Contruct the G matrix;
[nr, nc]=size(Zi);
G=2*eye(nr);
for i=1:1:nr-1;
    G(i,i+1)=-1;
    G(i+1,i)=-1;
end;

[rows, cols]=size(data_fd);
y=data_fd(:,cols);
x=[data_fd(:,1:4) dummy3]; %it does not matter whether dummy2 or dummy 3 is used. 
Z=[Z data_fd(:,2:4) dummy3]; %it does not matter whether dummy2 or dummy 3 is used. 

W=kron(eye(N),G);
Pz=Z*((Z'*W*Z)\Z');
beta_ab =(x'*Pz*x)\(x'*Pz*y);
[rowsy, colsy]=size(y);
[rowsx, colsx]=size(x); 


beta_ab_se=sum((y-x*beta_ah).^2)/(rowsy-colsx)*((x'*Pz*x)\eye(colsx));

%delta_y_hat=Pz*data_fd(:,1);
%x=[delta_y_hat data_fd(:,2:4) dummy2]
%beta_ab1=x\y;    %Check the 2SLS estimator


%% second step estimator

[nrow, ncol]=size(Z);
WW=zeros(ncol,ncol);
[~,TT]=size(G);
u_ab=y-x*beta_ab;
 for i=1:1:N;
     ind=(i-1)*TT+1:i*TT;
     Zi=Z(ind,:);
     WW=WW+Zi'*u_ab(ind)*(u_ab(ind))'*Zi;
 end; 
 
Pwz=Z*pinv(WW/N)*Z';

beta_ab_2step = (x'*Pwz*x)\(x'*Pz*y);



%% Arellano and Bond Estimator using lnc(t-2)as instruments 
Z=[];

T=30; 
for i=1:1:N;
 yi=lnc((i-1)*T+1:i*T,:);
 Zi=yi(1);
  for t=2:1:T-2;
    Zi=blkdiag(Zi,(yi(t))'); 
  end;
Z=[Z;Zi];
end;

y=data_fd(:,cols);
x=[data_fd(:,1:4) dummy3];

Z=[Z data_fd(:,2:4) dummy3];
Pz=Z*((Z'*W*Z)\(Z'));
beta_ab1= (x'*Pz*x)\(x'*Pz*y);


%% Arellano and Bond Estimator using lnc(t-2) lnc(t-3) as instruments 

Z=[];
Zi=[];

T=30; 
for i=1:1:N;
 yi=lnc((i-1)*T+1:i*T,:);
 Zi=yi(1);
  for t=2:1:T-2;
     Zi=blkdiag(Zi,(yi(t-1:t,1))'); 
  end;
Z=[Z;Zi];
end;

y=data_fd(:,cols);
x=[data_fd(:,1:4) dummy3];

Z=[Z data_fd(:,2:4) dummy3];
Pz=Z*((Z'*W*Z)\Z');
beta_ab2= (x'*Pz*x)\(x'*Pz*y);