%% Secrecy rate vs Power
clear; clc;
% seed = 100;
% rng(seed)

disp('--------------<Secrecy rate vs Power>--------------')

Power_dB = -5:2:15;
Power = 10.^(Power_dB/10);

% Number of Antennas
Nt = 8;  % no. of transmit antennas
Nr_1 = 1;  % no. of legitimate receive antennas
Nr_2 = 1;  % no. of eavesdropper receive antennas
M = 8;  % no. of IRS shifts

disp("No. of Antennas ")
disp("Tx: " + string(Nt))
disp("Rx: " + string(Nr_1))
disp("Eves: " + string(Nr_2))
disp("IRS: " + string(M))
disp('---------------------------------------------------')

psk_level = 8; % PSK order, bit = log2(psk_level);

Iter=3;

Rate_FS = zeros(length(Power_dB),1);
Rate_CD = zeros(length(Power_dB),1);
Rate_without_IRS = zeros(length(Power_dB),1);

for iter=1:Iter
    disp("current iter = "+string(iter))
    disp("current time: " + datestr(datetime('now')))
    disp('---------------------------------------------------')

    [H1, H2, R1, R2, T] = IRS_channel(Nt, Nr_1, Nr_2, M);
    
    for k=1:length(Power_dB)      
        % System without IRS
        [rate_no_IRS, ~, ~, ~] = MIMO_Capacity_with_fixed_PHI(H1, H2, Power(k));
        Rate_without_IRS(k) = Rate_without_IRS(k) + rate_no_IRS;
        
        % Proposed Algorithm
        [rate_CD, ~, ~, ~] = proposed_algorithm(H1, H2, R1, R2, T, Power(k));
        Rate_CD(k) = Rate_CD(k) + rate_CD;
        
        % Exhaustive search
        [rate_FS, ~, ~] = full_search(H1, H2, R1, R2, T, Power(k), psk_level);
        Rate_FS(k) = Rate_FS(k) + rate_FS;
    end
end
Rate_without_IRS = Rate_without_IRS/Iter;
Rate_CD = Rate_CD/Iter;
Rate_FS = Rate_FS/Iter;

figure;
plot(Power_dB, Rate_FS,'+-'); hold on
plot(Power_dB, Rate_CD,'s-'); 
plot(Power_dB, Rate_without_IRS,'*-');
legend('Exhaustive method','CD (Algorithm 1)', 'Without IRS','location','northwest'); 
xlabel('SNR (dBm)'); ylabel('Secrecy Rate (bps/Hz)'); grid on;

filename = "results\Secrecy rate vs Power _"+"Nt("+string(Nt)+")Nr_1("+string(Nr_1)+")Nr_2("+string(Nr_2)+")M("+string(M)+")Iter("+string(Iter)+")";
savefig(filename)
% save(filename)


%% Secrecy rate versus M
clear; clc;
% seed = 100;
% rng(seed)

disp('-----<Secrecy rate vs Number of IRS elements>------')

Power_dB = 10;
power = 10^(Power_dB/10);

% Number of Antennas
Nt = 8;  % no. of transmit antennas
Nr_1 = 1;  % no. of legitimate receive antennas
Nr_2 = 1;  % no. of eavesdropper receive antennas
M = [1, 4:4:32];  % no. of IRS shifts

disp("No. of Antennas ")
disp("Tx: " + string(Nt))
disp("Rx: " + string(Nr_1))
disp("Eves: " + string(Nr_2))
disp("Power: " + string(power) + " dB")
disp('---------------------------------------------------')

Iter=100;

Rate_without_IRS = 0;
Rate_CD = zeros(length(M),1);
time_without_IRS = 0;
time_CD = zeros(length(M),1);

for iter=1:Iter
    disp("current iter = "+string(iter))
    disp("current time: " + datestr(datetime('now')))
    disp('---------------------------------------------------')
    
    [H1, H2, R1, R2, T] = IRS_channel(Nt, Nr_1, Nr_2, max(M));
    
    tic % System without IRS
    [rate_no_IRS, ~, ~, ~] = MIMO_Capacity_with_fixed_PHI(H1, H2, power);
    time_without_IRS = time_without_IRS + toc;
    Rate_without_IRS = Rate_without_IRS + rate_no_IRS;
    
    for k=1:length(M)      
        tic % Proposed Algorithm
        [rate_CD, ~, ~, ~] = proposed_algorithm(H1, H2, R1(:,1:M(k)), R2(:,1:M(k)), T(:,1:M(k)), power);
        time_CD(k) = time_CD(k) + toc;
        Rate_CD(k) = Rate_CD(k) + rate_CD;
    end
end
Rate_without_IRS = Rate_without_IRS/Iter;
Rate_CD = Rate_CD/Iter;
time_without_IRS = time_without_IRS/Iter;
time_CD = time_CD/Iter;

% figure Secrecy rate versus M
figure;
plot(M, Rate_CD,'s-'); hold on
plot(M, Rate_without_IRS*ones(length(M),1),'*-')
legend('CD (Algorithm 1)','Without IRS','location','northwest')
xlabel('Number of IRS''s reflecting units'); ylabel('Secrecy Rate (bps/Hz)'); grid on

filename = "results\Secrecy rate vs M _"+"Nt("+string(Nt)+")Nr_1("+string(Nr_1)+")Nr_2("+string(Nr_2)+")P("+string(Power_dB)+")Iter("+string(Iter)+")";
savefig(filename)
% save(filename)

% figure time versus M
figure;
plot(M, time_CD,'s-'); hold on
plot(M, time_without_IRS*ones(length(M),1)','*-');
legend('CD (Algorithm 1)','Without IRS','location','northwest')
xlabel('Number of IRS''s reflecting units'); ylabel('Time (sec)'); grid on

filename = "results\Time vs M _"+"Nt("+string(Nt)+")Nr_1("+string(Nr_1)+")Nr_2("+string(Nr_2)+")P("+string(Power_dB)+")Iter("+string(Iter)+")";
savefig(filename)


%% Secrecy rate versus Eves location (dx)
clear; clc;
% seed = 100;
% rng(seed)

disp('-------<Secrecy rate versus Eves location>---------')

Power_dB = 10;
power = 10^(Power_dB/10);

% Number of Antennas
Nt = 8;  % no. of transmit antennas
Nr_1 = 1;  % no. of legitimate receive antennas
Nr_2 = 1;  % no. of eavesdropper receive antennas
M = 32;  % no. of IRS shifts

disp("No. of Antennas ")
disp("Tx: " + string(Nt))
disp("Rx: " + string(Nr_1))
disp("Eves: " + string(Nr_2))
disp("IRS: " + string(M))
disp('---------------------------------------------------')

dx = 20:20:140; % eavesdropper location

Iter=100;

Rate_without_IRS = zeros(length(dx),1);
Rate_CD = zeros(length(dx),1);

for iter=1:Iter
    disp("current iter = "+string(iter))
    disp("current time: " + datestr(datetime('now')))
    disp('---------------------------------------------------')

    for k=1:length(dx)
        
        [H1, H2, R1, R2, T] = IRS_channel(Nt, Nr_1, Nr_2, M, dx(k));
        
        % lower bound (system without IRS)
        [rate_no_IRS, ~, ~, ~] = MIMO_Capacity_with_fixed_PHI(H1, H2, power);
        Rate_without_IRS(k) = Rate_without_IRS(k) + rate_no_IRS;
               
        % Proposed Algorithm
        [rate_CD, ~, ~, ~] = proposed_algorithm(H1, H2, R1(:,1:M), R2(:,1:M), T(:,1:M), power);
        Rate_CD(k) = Rate_CD(k) + rate_CD;
    end
end
Rate_without_IRS = Rate_without_IRS/Iter;
Rate_CD = Rate_CD/Iter;

figure;
plot(dx, Rate_CD,'s-'); hold on
plot(dx, Rate_without_IRS,'*-');
legend('CD (Algorithm 1)','Without IRS','location','northwest')
grid on; axis([20 140 1 5.5]);
xlabel('Eavesdropper location'); ylabel('Secrecy Rate (bps/Hz)');

filename = "results\Secrecy rate vs dx _"+"Nt("+string(Nt)+")Nr_1("+string(Nr_1)+")Nr_2("+string(Nr_2)+")M("+string(M)+")Iter("+string(Iter)+")";
savefig(filename)
% save(filename)

