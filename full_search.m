function [Rate, Q, alpha] = full_search(H1, H2, R1, R2, T, Power, Mod)

Nr_1 = size(H1,1);  % no. of legitimate receive antennas
Nr_2 = size(H2,1);  % no. of eavesdropper receive antennas
M = size(R1, 2);  % no. of IRS shifts

I1=eye(Nr_1);
I2=eye(Nr_2);

set_digit = pskmod(0:(Mod)-1,(Mod),pi/Mod);

rate = zeros(1, Mod^M);

for iter=1:(Mod)^M

    a = dec2base(iter-1,Mod,M);
    for i = 1:M
        alpha(i,1) = set_digit(uint8(str2num(a(i)))+1);
    end
    
    H1_tilde = H1 + R1*diag(alpha)*T';
    H2_tilde = H2 + R2*diag(alpha)*T';
    
    [~, Q, ~, ~] = MIMO_Capacity_with_fixed_PHI(H1_tilde, H2_tilde, Power);

    rate(iter) = real(log2(det(I1+H1_tilde*Q*H1_tilde')) - log2(det(I2+H2_tilde*Q*H2_tilde')));
end

[Rate, index] = max(rate);

a = dec2base(index-1,Mod,M);
for i = 1:M
    alpha(i,1) = set_digit(uint8(str2num(a(i)))+1);
end



