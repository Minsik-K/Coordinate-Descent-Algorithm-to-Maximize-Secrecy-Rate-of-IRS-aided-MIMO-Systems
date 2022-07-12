function [rate, Q, alpha, rate_traj] = proposed_algorithm(H1, H2, R1, R2, T, Power)

Iter = 100;
eps=10^-3;

NT = size(H1,2);  % no. of transmit antennas
Nr_1 = size(H1,1);  % no. of legitimate receive antennas
Nr_2 = size(H2,1);  % no. of eavesdropper receive antennas
M = size(R1, 2);  % no. of IRS shifts

I1=eye(Nr_1);
I2=eye(Nr_2);

% Alpha initialization
[ ~, Q, ~, ~ ] = MIMO_Capacity_with_fixed_PHI( H1, H2, Power );
alpha = ones(M,1);
[~, alpha] = CD(H1, H2, R1, R2, T, Q, alpha, Iter, eps);

H1_tilde = H1 + R1*diag(alpha)*T';
H2_tilde = H2 + R2*diag(alpha)*T';

old_rate=0;
for iter=1:Iter

    [~, Q, ~, ~] = MIMO_Capacity_with_fixed_PHI(H1_tilde, H2_tilde, Power);
    
    [~, alpha] = CD(H1, H2, R1, R2, T, Q, alpha, Iter, eps);
        
    H1_tilde = H1 + R1*diag(alpha)*T';
    H2_tilde = H2 + R2*diag(alpha)*T';
    
    rate = real(log2(det(I1+H1_tilde*Q*H1_tilde')) - log2(det(I2+H2_tilde*Q*H2_tilde')));
               
    rate_traj(iter)=rate;
    if abs(rate-old_rate)<eps
        break
    else
        old_rate=rate;
    end    
end
[rate, Q, ~, ~] = MIMO_Capacity_with_fixed_PHI( H1_tilde, H2_tilde, Power); % last Q update


function [rate, alpha] = CD(H1, H2, R1, R2, T, Q, alpha, Iter, eps)

Nr_1 = size(H1,1);  % no. of legitimate receive antennas
Nr_2 = size(H2,1);  % no. of eavesdropper receive antennas
M = size(R1, 2);  % no. of IRS shifts

I1=eye(Nr_1);
I2=eye(Nr_2);

H1_ = H1 + R1*diag(alpha)*T';
H2_ = H2 + R2*diag(alpha)*T';

psk_level = 0; % PSK order, bit = log2(psk_level);

old_rate = 0;
for iter=1:Iter

    for m=1:M

        H1_ = H1_ - alpha(m)*R1(:,m)*T(:,m)';
        H2_ = H2_ - alpha(m)*R2(:,m)*T(:,m)';

        A1m = I1 + H1_*Q*H1_' + R1(:,m)*T(:,m)'*Q*T(:,m)*R1(:,m)';
        A2m = I2 + H2_*Q*H2_' + R2(:,m)*T(:,m)'*Q*T(:,m)*R2(:,m)';

        A1m_inv = inv(sqrtm(A1m));
        A2m_inv = inv(sqrtm(A2m));

        p1m = A1m_inv * R1(:,m)+10^-20;
        q1m = A1m_inv * H1_*Q*T(:,m)+10^-20;

        p2m = A2m_inv * R2(:,m)+10^-20;
        q2m = A2m_inv * H2_*Q*T(:,m)+10^-20;

        a1 = 2*(q1m'*p1m)' / (1 - (norm(p1m)^2*norm(q1m)^2 - abs(q1m'*p1m)^2));
        a2 = 2*(q2m'*p2m)' / (1 - (norm(p2m)^2*norm(q2m)^2 - abs(q2m'*p2m)^2));

        a=a1; b=a2;
        a_2=abs(a)^2; b_2=abs(b)^2;
        ab=real(a*b');
        lambda=(1-ab+sqrt( (1-ab)^2-(1-a_2)*(1-b_2))) / (1-b_2);
        alpha(m) = (a1-lambda*a2)/abs((a1-lambda*a2));
        
        if psk_level ~= 0  % Projection to PSK constellation
            alpha(m) = pskmod(pskdemod(alpha(m),psk_level,pi/psk_level),psk_level,pi/psk_level);
        end
    
        H1_ = H1_ + alpha(m)*R1(:,m)*T(:,m)';
        H2_ = H2_ + alpha(m)*R2(:,m)*T(:,m)';
    end    

    rate = real(log2(det(I1+H1_*Q*H1_')) - log2(det(I2+H2_*Q*H2_')));

    if rate-old_rate<eps
        break
    else
        old_rate=rate;
    end 
end