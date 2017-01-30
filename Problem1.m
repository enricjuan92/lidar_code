% PROBLEM 1: LINK-BUDGET (ELASTIC-BACKSCATTER LIDAR)
close all;
clear all;
clc;

%% PARAMETERS

    % LASER
    E = 160e-3; %Energy [J]
    
    % TELESCOPE
    D_tele = 0.2032;  % Celestron Schmidt-Cassegrain C-8 [m]
    dsh    = 0.06858; % Shade diameter [m]
    f      = 2;       % focal length [m]
 
    % CONSTANTS
    c = 2.99793e8; % speed light [m/s]
    
    % BACKGROUND-RADIANCE
    Lmoon =  3e-11; % Moon's radiance (full Moon) [W·cm-2·nm-1·sr-1]
    
    %INTERFERENCE FILTER
    Delta_Lambda = 10; % Bandwidth [nm]
    
    % PHOTODIODE
    D_apd = 3; % Active area diameter [mm]
    
    % ATMOSPHERE
    % Aerosol component
    VM        = 39.12; % Visibility margin (532 nm) [km]
    Alpha_aer = 3.912 / VM; 
    aer_ratio = 25; % Lidar ratio [sr]
    R_pbl     = 3; % Boundary-layer height, RPBL [km]
    % Molecular component (average)
    Alpha_mol = 0.01; % Rayleigh's extinction [km^-1]
    mol_ratio = 8 * pi / 3; % Rayleigh's ratio
    
%% 1. Determine the system constant, K(lambda) [W·km3]

Ar = pi * ((0.2032/2)^2-(0.06858/2)^2); % [m^2]
Ar_km = Ar * 1e-6; % [km^2]
c_km = c * 1e-3; % [km/s]
K = E * c_km * Ar_km / 2; % K = E*c*Ar/2 [W * km^3]
fprintf('[Answer Q1] K = %.2e W·Km^3\n', K); 

%% 2. Estimate the received background power under night-time operation, Pback

Ar_cm = Ar * 1e4; % [cm^2]
rd = D_apd * 1e-3 / 2; % [m]
FOV = atan(rd/f); % Field of View [rad]
Delta_Omega = pi*(sin(FOV))^2; % delta solid angle [sr]

Pback = Lmoon * Ar_cm * Delta_Omega * Delta_Lambda; % [W]

fprintf('[Answer Q2] Pback = %.2e W\n', Pback);

%% 3. Compute and plot the return power at the following ranges

Beta_aer = Alpha_aer / aer_ratio; % [km-1]
Beta_mol = Alpha_mol / mol_ratio; % [km-1]

R_km = [0.2 1 2 3 3.0001 4]; % range steps [km]
P_R  = zeros(1, length(R_km));

for i=1:length(R_km)
    P_R(i) = power_return(R_km(i), R_pbl, Alpha_aer, Alpha_mol, Beta_aer, Beta_mol, K); % [W]
end

% Plot
% figure
subplot(2, 2, 1);
semilogy(R_km,P_R,':r*');
title('Return Power P(R)')
xlabel('R [km]')
ylabel('P(R) [W]')
for i = 1:length(P_R)
    text(R_km(i), P_R(i), num2str(P_R(i)));
end

fprintf('[Answer Q3] Pback = [ %.2e %.2e %.2e %.2e %.2e %.2e ] W\n', P_R(1), P_R(2), P_R(3), P_R(4), P_R(5), P_R(6));

%% 4. Compute the receiver-chain voltage responsivity, Rv, and the net voltage responsivity Rv'

Gt = 5750; % [Ohms]
Gac = 20.3; % [V/V]
GT = Gt * Gac; % [Ohms]
Rio = 240e-3; % [A/W]
M = 150; % [·]
Rv = Rio * GT * M; % [V/W]
T1 = 0.6; % [·]
T2 = 0.65; % [·]
L = T1 * T2; % [·]
Rv_net = Rv * L; % [V/W]

fprintf('[Answer Q4] Rv = %.2e V/W and Rv_net = %.2e V/W\n', Rv, Rv_net);

%% 5. a) Compute the range-dependent signal-to-noise ratio (consider the ranges of question 3),SNR(R).
%%    b) Identify the noise-dominant system-operation mode.

q = 1.602e-19; % [C]
F = 4.5; % [·]
Ids = 7.64e-8; % [A]
Idb = 3.10e-10; % [A]
Sigma_thi = 5e-12; % [A/sqrt(Hz)]
B = 10e6; % [Hz]

Sigma_shs = zeros(1, length(R_km));
Sigma_v   = zeros(1, length(R_km));
Sigma_shd = zeros(1, length(R_km));
Sigma_th  = zeros(1, length(R_km));
N         = zeros(1, length(R_km));
SNR       = zeros(1, length(R_km));
SNR_db    = zeros(1, length(R_km));

for i=1:length(R_km)
    Sigma_shs(i) = 2 * q * GT^2 * F * M^2 * Rio * (P_R(i) + Pback) * L; % [V^2/Hz]
    Sigma_shd(i) = 2 * q * GT^2 * (Ids + F * M^2 * Idb); % [V^2/Hz]
    Sigma_th(i)  = Sigma_thi^2 * GT^2; % [V^2/Hz]
    Sigma_v(i)   = Sigma_shs(i) + Sigma_shd(i) + Sigma_th(i);
    N(i)         = Sigma_v(i) * B;
    SNR(i)       = Rv * L * P_R(i) / sqrt(N(i));
    SNR_db(i)    = 20 * log10(SNR(i));
end

% Plot SNR
% figure
subplot(2, 2, 2);
plot(R_km, SNR_db,':r*');
title('Signal to Noise Ratio SNR(R) [dB]')
xlabel('R [km]')
ylabel('SNR(R) [dB]') 
for i = 1:length(SNR_db)
    text(R_km(i), SNR_db(i), num2str(SNR_db(i)));
end

% Plot Noises
% figure
subplot(2, 2, 3);
semilogy(R_km, Sigma_shs * B,':r*');
title('N(V^2)')
xlabel('R [km]')
ylabel('N(V^2)') 
for i = 1:length(Sigma_shs)
    text(R_km(i), Sigma_shs(i) * B, num2str(Sigma_shs(i) * B));
end
hold on
semilogy(R_km, Sigma_shd * B,':b*');
text(R_km(1), Sigma_shd(1) * B, num2str(Sigma_shd(1) * B));
hold on
semilogy(R_km, Sigma_th * B,':g*');
text(R_km(1), Sigma_th(1) * B, num2str(Sigma_th(1) * B));
legend('\sigma^2_{shs}','\sigma^2_{shd}','\sigma^2_{th}')

%% 6. Assess the approximate laser-radar maximum range (SNR(Rmax)=1)

R_kmQ6 = 0.2:0.01:20;

P_R_Q6       = zeros(1, length(R_kmQ6));
Sigma_shs_Q6 = zeros(1, length(R_kmQ6));
Sigma_v_Q6   = zeros(1, length(R_kmQ6));
Sigma_shd_Q6 = zeros(1, length(R_kmQ6));
Sigma_th_Q6  = zeros(1, length(R_kmQ6));
N_Q6         = zeros(1, length(R_kmQ6));
SNR_Q6       = zeros(1, length(R_kmQ6));
SNR_db_Q6    = zeros(1, length(R_kmQ6));

for i=1:length(R_kmQ6)
   P_R_Q6(i)       = power_return(R_kmQ6(i), R_pbl, Alpha_aer, Alpha_mol, Beta_aer, Beta_mol, K);
   Sigma_shs_Q6(i) = 2 * q * GT^2 * F * M^2 * Rio * (P_R_Q6(i) + Pback) * L; % [V^2/Hz]
   Sigma_v_Q6(i)   = Sigma_shs_Q6(i) + Sigma_shd(1) + Sigma_th(1);
   N_Q6(i)         = Sigma_v_Q6(i) * B; % [V^2]
   SNR_Q6(i)       = (Rv * L * P_R_Q6(i)) / (sqrt(N_Q6(i)));
   SNR_db_Q6(i)    = 20 * log10(SNR_Q6(i));
end

R_max = R_kmQ6(find(SNR_Q6 >= 1, 1, 'last' ));

fprintf('[Answer Q6] Rmax = %.3f Km\n', R_max);

%% 7. How many pulses are needed to integrate in order to ensure a SNRv (voltage signal-to-noise ratio) of 40 dB at 3-km range? What is the resulting observation time of the lidar instrument?

SNR_40dB = 10^(40/20);
Npulses = ceil(SNR_40dB^2 / SNR(4)^2);
freq = 20; % [Hz]
Tobs = (1/20) * Npulses;

fprintf('[Answer Q7] Npulses = %d \n', Npulses);

%% 8. Now, consider a Raman system of similar specs. If for Raman systems the return signal is typically 3 orders of magnitude lower than for their elastic system counterparts, discuss on the feasibility of day-time operation.

R_raman = linspace(0.2, 4);
P_raman = 1e-3*(K./R_raman.^2).*(Beta_aer+Beta_mol).*exp(-2*(Alpha_aer+Alpha_mol).*R_raman); % [W]

% figure
subplot(2, 2, 4);
semilogy(R_raman,P_raman,'-');
title('Return Raman Power Praman(R)')
xlabel('R [km]')
ylabel('Praman(R) [W]') 
hold on
semilogy(R_km, P_R, ':r*');
hold on
semilogy(R_raman,ones(size(R_raman))*3e-6,'g');
legend('Raman','Elastic','Day-time background')

%% 9. Compute the photodiode NEP and its quantum efficiency (eta)

h = 6.6262e-34;
Lambda = 532e-9;

NEP = (sqrt(Sigma_shd(1) / GT^2))/(Rio * M); % [W/sqrt(Hz)]
Eta = (Rio * h * c)/(q * Lambda);

fprintf('[Answer Q9] NEP = %.2e W/sqrt(Hz) and Eta = %.2f\n', NEP, Eta*100);


% 10 Compute the system NEP (NEPs)

NEPs = (sqrt(Sigma_shd(1) + Sigma_th(1)))/(Rv_net); % [W/sqrt(Hz)]

fprintf('[Answer Q10] NEP = %.2e W/sqrt(Hz)\n', NEPs);
