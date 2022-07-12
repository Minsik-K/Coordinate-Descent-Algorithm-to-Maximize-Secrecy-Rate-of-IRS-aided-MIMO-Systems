function [H1, H2, R1, R2, T] = IRS_channel(Nt, N1, N2, M, d_x)
%
% Our channel model is based on following paper.
% M. Cui, G. Zhang, and R. Zhang, ``Secure wireless communication via intelligent reflecting 
% surface," \emph{IEEE Wireless Commun. Lett.}, vol. 8, no. 5, pp. 1410-1414, Oct. 2019.
%

option=0;   % Rician fiading LOS components

if nargin<5
    d_x=140;   % Rician fiading LOS components
end   

zeta0 = 10^-3;
d0 = 1;

% The path loss exponent
beta_TR = 3.6;
beta_TI = 2.2;
beta_TE = 2.2;
beta_IR = 2.2;
beta_IE = 2.2;

alpha_AU=3.0; alpha_AE=3.0;
alpha_IU=1.5; alpha_IE=1.5; alpha_AI=1.5;

tx_location = [0, 0];
% rx_location = [150, 0];
% eves_location = [d_x, 0];
rx_location = [150+randn(1,1), 0];
eves_location = [d_x+randn(1,1), 0];
irs_location = [10, 10];

d_AU = norm(tx_location-rx_location);   % d_TR
d_AI = norm(tx_location-irs_location);  % d_TI
d_AE = norm(tx_location-eves_location);  % d_TE

d_IU = norm(irs_location-rx_location);  % d_IR
d_IE = norm(irs_location-eves_location);

% d_AU=150;
% d_AE=145;
% d_IE=5;
% d_IU=sqrt((d_x)^2+d_IE^2);
% d_AI=sqrt((d_AU-d_x)^2+d_IE^2);

sigma=sqrt(10^-8);  % noise variance -80 dBm

R_c=zeros(Nt,Nt);
for i=1:Nt
    for j=1:Nt
        R_c(i,j)=0.95^(abs(i-j));
    end
end 

gamma = 1;
R1_review=zeros(Nt,Nt);
% D = lambda/2
for i=1:Nt
    for j=1:Nt
        u_n = [0, mod(i-1, Nt/4), floor((i-1) / 4)];
        u_m = [0, mod(j-1, Nt/4), floor((j-1) / 4)];
%         R_review(i,j) = gamma * sin( 2*pi * norm(u_n-u_m) * D/lambda ) / (2*pi * norm(u_n-u_m) * D/lambda);
        R1_review(i,j) = gamma * sinc(norm(u_n-u_m));

    end
end 

R2_review=zeros(M,M);
% D = lambda/2
for i=1:M
    for j=1:M
        u_n = [0, mod(i-1, M/4), floor((i-1) / 4)];
        u_m = [0, mod(j-1, M/4), floor((j-1) / 4)];
%         R_review(i,j) = gamma * sin( 2*pi * norm(u_n-u_m) * D/lambda ) / (2*pi * norm(u_n-u_m) * D/lambda);
        R2_review(i,j) = gamma * sinc(norm(u_n-u_m));

    end
end 

K=1;
LOS = sqrt(K/(K+1));
NLOS = sqrt(1/(K+1));

%c_L = randn(1,Nt)+1i*randn(1,Nt); c_L=c_L./abs(c_L);
if option==0
    c_L = ones(N1,Nt);               % Rician fiading LOS components with same phase
else    
	c_L = exp(1i*rand(N1,Nt)*2*pi);  % Rician fiading LOS components with random phase
end    
c_NL = (randn(N1,Nt)+1i*randn(N1,Nt))/sqrt(2);
h_AU = sqrt(zeta0*(d0/d_AU)^alpha_AU) * (LOS*c_L + NLOS*c_NL*sqrtm(R1_review)) / sigma;

%c_L = randn(1,Nt)+1i*randn(1,Nt); c_L=c_L./abs(c_L);
if option==0
    c_L = ones(N2,Nt);
else    
	c_L = exp(1i*rand(N2,Nt)*2*pi);
end    
c_NL = (randn(N2,Nt)+1i*randn(N2,Nt))/sqrt(2);
h_AE = sqrt(zeta0*(d0/d_AE)^alpha_AE) * (LOS*c_L + NLOS*c_NL*sqrtm(R1_review)) / sigma;

%c_L = randn(1,M)+1i*randn(1,M); c_L=c_L./abs(c_L);
if option==0
    c_L = ones(N1,M);
else    
	c_L = exp(1i*rand(N1,M)*2*pi);
end    
c_NL = (randn(N1,M)+1i*randn(N1,M))/sqrt(2);
h_IU = sqrt(zeta0*(d0/d_IU)^alpha_IU) * (LOS*c_L + NLOS*c_NL*sqrtm(R2_review)) / sigma;

%c_L = randn(1,M)+1i*randn(1,M); c_L=c_L./abs(c_L);
if option==0
    c_L = ones(N2,M);
else    
	c_L = exp(1i*rand(N2,M)*2*pi);
end    
c_NL = (randn(N2,M)+1i*randn(N2,M))/sqrt(2);
h_IE = sqrt(zeta0*(d0/d_IE)^alpha_IE) * (LOS*c_L + NLOS*c_NL*sqrtm(R2_review)) / sigma;

%c_L = randn(M,Nt)+1i*randn(M,Nt); c_L=c_L./abs(c_L);
if option==0
    c_L = ones(M,Nt);
else    
	c_L = exp(1i*rand(M,Nt)*2*pi);
end    
c_NL = (randn(M,Nt)+1i*randn(M,Nt))/sqrt(2);
H_AI = sqrt(zeta0*(d0/d_AI)^alpha_AI) * (LOS*c_L + NLOS*c_NL*sqrtm(R1_review));  % Division by sigma is not necessary.

H1=h_AU; H2=h_AE;
R1=h_IU; R2=h_IE;
T=H_AI';
end

