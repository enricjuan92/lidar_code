%Problem 2
% close all;
clear all;
clc;

%% PARAMETERS


%% 1. Determine the system constant, K [W·km3]

c = 2.99793e8; % [m/s]
E = 150e-3; % [J]
d_primary = 0.3; % [m]
R = (d_primary/2)^2-(0.06858/2)^2; % [m]
Ar = pi*R; % [m^2]
K_1 = ((E*c)/2)*Ar; % [Wm^3]
K = K_1*1e-9; % [Wkm^3]

fprintf('[Answer Q1] K = %.2e W·Km^3\n', K);

%% 2. Estimate the received background power for both the elastic and Raman channels (Pback,0 and Pback,R, respectively) under (see “day-time/night-time” parameter) operation.

Lsun = 3e-6; % [Wcm-2nm-1sr-1]
Ar_1 = Ar*1e4; % [cm^2]
rd = 3e-3/2; % [m]
f = 2; % [m]
FOV = atan(rd/f); % [rad]
Delta_Omega = pi*(sin(FOV))^2; % [sr]
Delta_Lambda_elastic = 10; %[nm]
Delta_Lambda_raman = 3; %[nm]
Pback_elastic = Lsun*Ar_1*Delta_Omega*Delta_Lambda_elastic; % [W]
Pback_raman = Lsun*Ar_1*Delta_Omega*Delta_Lambda_raman; % [W]

fprintf('[Answer Q2] Pback(elastic) = %.2e W and Pback(Raman) = %.2e W\n', Pback_elastic, Pback_raman);

%% 3. Plot both elastic and Raman-return powers (P0(R) and PR(R), respectively). Superimpose Pback,0 and Pback,R plots from question 2 results.

R_15km      = 0:0.1:15;
R_PBL       = 3; % [km]]
VM          = 39.12;
Alpha_aer   = 3.912/VM; % [km-1]
Beta_aer    = Alpha_aer/25; % [km-1]

for i=1:151
   Alpha_mol = 1.2569e-2 - 7.7599e-4*R_15km(i); % [km-1]
   Beta_mol  = Alpha_mol/(8*pi/3); % [km-1]
   p_r(i) = power_return(R_15km(i), R_PBL, Alpha_aer, Alpha_mol, Beta_aer, Beta_mol, K);
   NR = 2.1145e34 - 2.0022e33*R_15km(i) + 5.4585e31*R_15km(i)^2;
   RCross_section = 3.71e-41; % [km2sr-1]
   Lambda_r  = 607.4e-9; % [m]
   Lambda_0  = 532e-9; % [m]
   Alpha_mol_int = 1.2569e-2*R_15km(i) - 7.7599e-4*R_15km(i)*R_15km(i)/2;
   Alpha_raman_aer = ((Lambda_r/Lambda_0)^(-1.8))*Alpha_aer;
   Alpha_raman_mol = 7.3219e-3 - 4.5204e-4*R_15km(i); % [km-1]
   Alpha_raman_mol_int = 7.3219e-3*R_15km(i) - 4.5204e-4*R_15km(i)*R_15km(i)/2;
   p_r_raman(i) = (K/R_15km(i))*(NR*RCross_section)*exp(-2 * (Alpha_aer*R_15km(i) + Alpha_mol_int*+Alpha_raman_aer* R_15km(i)+Alpha_raman_mol_int) );
end

figure
subplot(2,2,1);
semilogy(R_15km,p_r,'r');
for i = 1:length(p_r)
    if R_15km(i) == 0.2 || R_15km(i) == 3  || R_15km(i) == 8
        label = strcat('\leftarrow  X: ', num2str(R_15km(i)),' km', ' Y: ', num2str(p_r(i)), ' W');
        text(R_15km(i), p_r(i), label);
    end
end
hold on
semilogy(R_15km,p_r_raman,'b');
for i = 1:length(p_r_raman)
    if R_15km(i) == 3  || R_15km(i) == 8
        label = strcat('\leftarrow  X: ', num2str(R_15km(i)),' km', ' Y: ', num2str(p_r_raman(i)), ' W');
        text(R_15km(i), p_r_raman(i), label);
    end
end
hold on
semilogy(R_15km,ones(size(R_15km))*Pback_elastic,'c');
hold on
semilogy(R_15km,ones(size(R_15km))*Pback_raman,'g');
title('Return Power P(R)')
xlabel('R [km]')
ylabel('P(R) [W]') 
legend('Elastic','Raman','Pback Elastic','Pback Raman')


%% 4. Compute, for the elastic and Raman channels, receiver-chain voltage responsivities (Rv,0 and Rv,R, respectively), and net voltage responsivities (i.e., including spectral optical losses; Rv,0’ and Rv,R’, respectively).

Gt              = 5750; % [Ohms]
Gac             = 20.3; % [V/V]
GT_elastic      = Gt * Gac; % [Ohms]
Rio_elastic     = 240e-3; % [A/W]
M_elastic       = 150; % [·]
Rv_elastic      = Rio_elastic * GT_elastic * M_elastic; % [V/W]
T1              = 0.6; % [·]
T2              = 0.65; % [·]
L               = T1 * T2; % [·]
Rv_net_elastic  = Rv_elastic * L; % [V/W]
M_raman         = 3e6; % [·]
Ri              = 3e4; % [A/W]
Rio_raman       = Ri/M_raman; % [A/W]
GT_raman        = 50; % [Ohms] 
Rv_raman        = Rio_raman*GT_raman*M_raman; % [V/W]
Rv_net_raman    = Rv_raman*L; % [V/W]

fprintf('[Answer Q4] Rv(elastic) = %.2e V/W and Rv_net(elastic) = %.2e V/W\n', Rv_elastic, Rv_net_elastic);
fprintf('[Answer Q4] Rv(Raman) = %.2e V/W and Rv_net(Raman) = %.2e V/W\n', Rv_raman, Rv_net_raman);

%% 5. a) Assuming analog detection, plot the elastic range-dependent signal-to-noise ratio, SNR0(R), 
%% at the output of the receiver chain (i.e., voltage ratio) and related shot photo-induced, shot-dark and thermal variances

p_r       = zeros(1, length(R_15km));
Sigma_shs = zeros(1, length(R_15km));
Sigma_v   = zeros(1, length(R_15km));
N         = zeros(1, length(R_15km));
SNR       = zeros(1, length(R_15km));
SNR_db    = zeros(1, length(R_15km));

R_15km      = 0:0.1:15;
q           = 1.602e-19; % [C]
F           = 4.5; % [·]
Ids         = 7.64e-8; % [A]
Idb         = 3.10e-10; % [A]
Sigma_thi   = 5e-12; % [A/sqrt(Hz)]
B           = 10e6; % [Hz]

Alpha_aer = 1; % [km-1]
Beta_aer  = Alpha_aer/25; % [km-1]

for i=1:length(R_15km)
   Alpha_mol = 1.2569e-2 - 7.7599e-4 * R_15km(i); % [km-1]
   Beta_mol  = Alpha_mol/(8*pi/3); % [km-1]
   
   p_r(i) = power_return(R_15km(i), R_PBL, Alpha_aer, Alpha_mol, Beta_aer, Beta_mol, K);
   
   Sigma_shs(i) = 2 * q * GT_elastic^2 * F * M_elastic^2 * Rio_elastic * (p_r(i) + Pback_elastic) * L; % [V^2/Hz]
   Sigma_shd = 2 * q * GT_elastic^2 * (Ids + F * M_elastic^2 * Idb); % [V^2/Hz]
   Sigma_th = Sigma_thi^2 * GT_elastic^2; % [V^2/Hz]
   Sigma_v(i) = Sigma_shs(i) + Sigma_shd + Sigma_th;
  
   N(i) = Sigma_v(i) * B; % [V^2]
   SNR(i) = (Rv_elastic * L * p_r(i))/(sqrt(N(i)));
   SNR_db(i) = 20 * log10(SNR(i));
end

% Plot Noises
subplot(2, 2, 2);
semilogy(R_15km, Sigma_shs, 'r');
title('N(V^2)')
xlabel('R [km]')
ylabel('N(V^2)') 
hold on
semilogy(R_15km, ones(size(R_15km)) * Sigma_shd, 'b');
hold on
semilogy(R_15km, ones(size(R_15km)) * Sigma_th, 'g');
legend('\sigma^2_{shs}','\sigma^2_{shd}','\sigma^2_{th}')

% Plot SNR
subplot(2,2,3);
plot(R_15km, SNR_db,'r');
for i = 1:length(SNR_db)
    if R_15km(i) == 0.2 || R_15km(i) == 1 || R_15km(i) == 2 ...
    || R_15km(i) == 3   || R_15km(i) == 4 || R_15km(i) == 8 ...
    || R_15km(i) == 13
        label = strcat('\leftarrow  X: ', num2str(R_15km(i)),' km', ' Y: ', num2str(SNR_db(i)), ' dB');
        text(R_15km(i), SNR_db(i), label);
    end
end

title('Signal to Noise Ratio SNR(R) [dB]')
xlabel('R [km]')
ylabel('SNR(R) [dB]') 
