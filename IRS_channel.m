function [H1, H2, R1, R2, T] = IRS_channel(Nt, N1, N2, M, d_x)
%
% Our channel model is based on following paper.
% M. Cui, G. Zhang, and R. Zhang, ``Secure wireless communication via intelligent reflecting 
% surface," \emph{IEEE Wireless Commun. Lett.}, vol. 8, no. 5, pp. 1410-1414, Oct. 2019.
%

option=0;   % Rician fiading LOS components

if nargin<5
    d_x=140;   % IRS location
end   

zeta0 = 10^-3;
d0 = 1;

% The path loss exponent
beta_TU=3.0; beta_TE=3.0;
beta_IU=1.5; beta_IE=1.5; 
beta_TI=1.5;

% location
tx_location = [0, 0];
rx_location = [150+randn(1,1), 0];
eves_location = [d_x+randn(1,1), 0];
irs_location = [10, 10];

% Distance
d_TR = norm(tx_location-rx_location);
d_TI = norm(tx_location-irs_location);
d_TE = norm(tx_location-eves_location);
d_IR = norm(irs_location-rx_location);
d_IE = norm(irs_location-eves_location);


sigma=sqrt(10^-8);  % noise variance -80 dBm

% Rayleigh fading model given in [17]
gamma = 1;
R1_review=zeros(Nt,Nt);
for i=1:Nt
    for j=1:Nt
        u_n = [0, mod(i-1, Nt/4), floor((i-1) / 4)];
        u_m = [0, mod(j-1, Nt/4), floor((j-1) / 4)];
        R1_review(i,j) = gamma * sinc(norm(u_n-u_m));
    end
end 

R2_review=zeros(M,M);
for i=1:M
    for j=1:M
        u_n = [0, mod(i-1, M/4), floor((i-1) / 4)];
        u_m = [0, mod(j-1, M/4), floor((j-1) / 4)];
        R2_review(i,j) = gamma * sinc(norm(u_n-u_m));
    end
end 

K=1; % Rician factor
LOS = sqrt(K/(K+1));
NLOS = sqrt(1/(K+1));

% -------------------- Tx - xxx channels ------------------------
% -------------------- Tx - USER channels -----------------------
if option==0
    c_L = ones(N1,Nt);               % Rician fiading LOS components with same phase
else    
	c_L = exp(1i*rand(N1,Nt)*2*pi);  % Rician fiading LOS components with random phase
end
c_NL = (randn(N1,Nt)+1i*randn(N1,Nt))/sqrt(2);
h_TU = sqrt(zeta0*(d0/d_TR)^beta_TU) * (LOS*c_L + NLOS*c_NL*sqrtm(R1_review)) / sigma;

% -------------------- Tx - IRS channels -----------------------
if option==0
    c_L = ones(M,Nt);
else    
	c_L = exp(1i*rand(M,Nt)*2*pi);
end    
c_NL = (randn(M,Nt)+1i*randn(M,Nt))/sqrt(2);
H_TI = sqrt(zeta0*(d0/d_TI)^beta_TI) * (LOS*c_L + NLOS*c_NL*sqrtm(R1_review));

% -------------------- Tx - Eves channels -----------------------
if option==0
    c_L = ones(N2,Nt);
else    
	c_L = exp(1i*rand(N2,Nt)*2*pi);
end    
c_NL = (randn(N2,Nt)+1i*randn(N2,Nt))/sqrt(2);
h_TE = sqrt(zeta0*(d0/d_TE)^beta_TE) * (LOS*c_L + NLOS*c_NL*sqrtm(R1_review)) / sigma;


% -------------------- IRS - xxx channels -----------------------
% ------------------- IRS - USER channels -----------------------
if option==0
    c_L = ones(N1,M);
else    
	c_L = exp(1i*rand(N1,M)*2*pi);
end    
c_NL = (randn(N1,M)+1i*randn(N1,M))/sqrt(2);
h_IU = sqrt(zeta0*(d0/d_IR)^beta_IU) * (LOS*c_L + NLOS*c_NL*sqrtm(R2_review)) / sigma;

% ------------------- IRS - USER channels -----------------------
if option==0
    c_L = ones(N2,M);
else    
	c_L = exp(1i*rand(N2,M)*2*pi);
end    
c_NL = (randn(N2,M)+1i*randn(N2,M))/sqrt(2);
h_IE = sqrt(zeta0*(d0/d_IE)^beta_IE) * (LOS*c_L + NLOS*c_NL*sqrtm(R2_review)) / sigma;

% ---------------------------------------------------------------


H1=h_TU; H2=h_TE;
R1=h_IU; R2=h_IE;
T=H_TI';
end

