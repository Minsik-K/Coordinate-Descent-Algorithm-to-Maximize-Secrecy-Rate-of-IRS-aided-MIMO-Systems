function [rate, S1, A1, e1 ] = MIMO_Capacity_with_fixed_PHI(H1, H2, Power)
%
% This algorithm is based on following paper.
% Q. Li, M. Hong, H.-T. Wai, Y.-F. Liu, W.-K. Ma, and Z.-Q. Luo, ``Transmit 
% solutions for MIMO wiretap channels using alternating optimization,?? 
% \emph{IEEE J. Sel. Areas Commun.}, vol. 31, no. 9, pp. 1714--1727, Sep. 2013.
%

M=size(H1,2);
Nr_1=size(H1,1);
Nr_2=size(H2,1); 

I1=eye(Nr_1);
I2=eye(Nr_2);

e1 = 0;

if Nr_1==1 && Nr_2==1
    A = eye(M)+Power*H1'*H1;
    B = eye(M)+Power*H2'*H2;
    [V, D]=eig(inv(B)*A);
    [~, index]=max(diag(D));
    e1=V(:,index);
    S1 = Power*e1*e1';
    A1 = 0; lambda = 0;
    rate = real(log2(det(I1+H1*S1*H1')) - log2(det(I2+H2*S1*H2')));
    return;
end    

S1=zeros(M,M); 
S2=zeros(M,M); 

S1=randn(M,M)+sqrt(-1)*randn(M,M);
S1=S1*S1';
S1=Power*S1/trace(S1);

epsilon1=10^-12;
epsilon2=10^-4;
gamma=1;

lambda_min=0;
lambda_max=5;

while(1)
    lambda=(lambda_min+lambda_max)/2;
    rate=0;

    while(1)
        rate_old=rate;
        if rcond(I2+H2*S1*H2')<10^-15
            disp('low rcond');
        end
        
        A1=H2'*inv(I2+H2*S1*H2')*H2;
        A1=(A1+A1')/2;

        S1_ = MIMO_Capacity_geneig(1, H1, I1, A1, lambda);
        S1 = (1-gamma)*S1 + gamma*S1_;
        S1=(S1+S1')/2;

        rate = real(log2(det(I1+H1*S1*H1')) - log2(det(I2+H2*S1*H2')));
        if abs(rate-rate_old)<epsilon2
            break;
        end
    end

    if real(trace(S1))>Power
        lambda_min=lambda;
    else
        lambda_max=lambda;
    end

    if abs(lambda_max-lambda_min) < epsilon1
        break;
    end
end

rate=real(rate);
end