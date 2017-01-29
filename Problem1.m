%Problem 1
close all;
clear all;
% 1. Determine the system constant, K(lambda) [W·km3]

c = 2.99793e8; % [m/s]
E = 160e-3; % [J]
R = (0.2032/2)^2-(0.06858/2)^2; % [m]
Ar = pi*R; % [m^2]
K_1 = ((E*c)/2)*Ar; % [Wm^3]
K = K_1*1e-9 % [Wkm^3]

% 2. Estimate the received background power under night-time operation, Pback

Lmoon = 3e-11; % [Wcm-2nm-1sr-1]
Ar_1 = Ar*1e4; % [cm^2]
rd = 3e-3/2; % [m]
f = 2; % [m]
FOV = atan(rd/f); % [rad]
Delta_Omega = pi*(sin(FOV))^2; % [sr]
Delta_Lambda = 10; %[nm]
Pback = Lmoon*Ar_1*Delta_Omega*Delta_Lambda % [W]

% 3. Compute and plot the return power at the following ranges

Alpha_aer = 0.1; % [km-1]
Alpha_mol = 0.01; % [km-1]
Beta_aer = Alpha_aer/25; % [km-1]
Beta_mol = Alpha_mol/(8*pi/3); % [km-1]

% 0.2 Km

R_200 = 0.2; % [km]
P_200 = (K/R_200^2)*(Beta_aer+Beta_mol)*exp(-2*(Alpha_aer+Alpha_mol)*R_200);

% 1 Km

R_1 = 1; % [km]
P_1= (K/R_1^2)*(Beta_aer+Beta_mol)*exp(-2*(Alpha_aer+Alpha_mol)*R_1); % [W]

% 2 Km

R_2 = 2; % [km]
P_2= (K/R_2^2)*(Beta_aer+Beta_mol)*exp(-2*(Alpha_aer+Alpha_mol)*R_2); % [W]

% 3 Km (-)

R_3d = 3; % [km]
P_3d= (K/R_3d^2)*(Beta_aer+Beta_mol)*exp(-2*(Alpha_aer+Alpha_mol)*R_3d); % [W]

% 3 Km (+)

R_3u = 3; % [km]
R_PBL = 3; % [km]
P_3u= (K/R_3u^2)*(Beta_mol)*exp(-2*(Alpha_aer+Alpha_mol)*R_PBL)*exp(-2*Alpha_mol*(R_3u-R_PBL)); % [W]

% 4 Km

R_4 = 4; % [km]
R_PBL = 3; % [km]
P_4= (K/R_4^2)*(Beta_mol)*exp(-2*(Alpha_aer+Alpha_mol)*R_PBL)*exp(-2*Alpha_mol*(R_4-R_PBL)); % [W]

% Plot

Rvector = [R_200;R_1;R_2;R_3d;R_3u;R_4];
Pvector = [P_200;P_1;P_2;P_3d;P_3u;P_4];
figure
semilogy(Rvector,Pvector,':r*');
title('Return Power P(R)')
xlabel('R [km]')
ylabel('P(R) [W]') 

% 4. Compute the receiver-chain voltage responsivity, Rv, and the net voltage responsivity Rv'

Gt = 5750; % [Ohms]
Gac = 20.3; % [V/V]
GT = Gt*Gac; % [Ohms]
Rio = 240e-3; % [A/W]
M = 150; % [·]
Rv = Rio*GT*M % [V/W]
T1 = 0.6; % [·]
T2 = 0.65; % [·]
L = T1 * T2; % [·]
Rv_net = Rv*L % [V/W]

% 5. a) Compute the range-dependent signal-to-noise ratio (consider the ranges of question 3),SNR(R).
%    b) Identify the noise-dominant system-operation mode.

% 0.2 Km
q = 1.602e-19; % [C]
F = 4.5; % [·]
Ids = 7.64e-8; % [A]
Idb = 3.10e-10; % [A]
Sigma_thi = 5e-12; % [A/sqrt(Hz)]
Sigma_shs_200 = 2*q*GT^2*F*M^2*Rio*(P_200+Pback)*L; % [V^2/Hz]
Sigma_shd = 2*q*GT^2*(Ids + F*M^2*Idb); % [V^2/Hz]
Sigma_th = Sigma_thi^2*GT^2; % [V^2/Hz]
Sigma_v_200 = Sigma_shs_200+Sigma_shd+Sigma_th;
B = 10e6; % [Hz]
N_200 = Sigma_v_200*B; % [V^2]
SNR_200 = (Rv*L*P_200)/(sqrt(N_200))
SNR_200_db = 20*log10(SNR_200)

% 1 Km
Sigma_shs_1 = 2*q*GT^2*F*M^2*Rio*(P_1+Pback)*L; % [V^2/Hz]
Sigma_v_1 = Sigma_shs_1+Sigma_shd+Sigma_th;
N_1 = Sigma_v_1*B; % [V^2]
SNR_1 = (Rv*L*P_1)/(sqrt(N_1))
SNR_1_db = 20*log10(SNR_1)

% 2 Km
Sigma_shs_2 = 2*q*GT^2*F*M^2*Rio*(P_2+Pback)*L; % [V^2/Hz]
Sigma_v_2 = Sigma_shs_2+Sigma_shd+Sigma_th;
N_2 = Sigma_v_2*B; % [V^2]
SNR_2 = (Rv*L*P_2)/(sqrt(N_2))
SNR_2_db = 20*log10(SNR_2)

% 3 Km (-)
Sigma_shs_3d = 2*q*GT^2*F*M^2*Rio*(P_3d+Pback)*L; % [V^2/Hz]
Sigma_v_3d = Sigma_shs_3d+Sigma_shd+Sigma_th;
N_3d = Sigma_v_3d*B; % [V^2]
SNR_3d = (Rv*L*P_3d)/(sqrt(N_3d))
SNR_3d_db = 20*log10(SNR_3d)

% 3 Km (+)
Sigma_shs_3u = 2*q*GT^2*F*M^2*Rio*(P_3u+Pback)*L; % [V^2/Hz]
Sigma_v_3u = Sigma_shs_3u+Sigma_shd+Sigma_th;
N_3u = Sigma_v_3u*B; % [V^2]
SNR_3u = (Rv*L*P_3u)/(sqrt(N_3u))
SNR_3u_db = 20*log10(SNR_3u)

% 4 Km (+)
Sigma_shs_4 = 2*q*GT^2*F*M^2*Rio*(P_4+Pback)*L; % [V^2/Hz]
Sigma_v_4 = Sigma_shs_4+Sigma_shd+Sigma_th;
N_4 = Sigma_v_4*B; % [V^2]
SNR_4 = (Rv*L*P_4)/(sqrt(N_4))
SNR_4_db = 20*log10(SNR_4)
% Plot SNR

Rvector = [R_200;R_1;R_2;R_3d;R_3u;R_4];
SNRvector = [SNR_200_db;SNR_1_db;SNR_2_db;SNR_3d_db;SNR_3u_db;SNR_4_db];
figure
plot(Rvector,SNRvector,':r*');
title('Signal to Noise Ratio SNR(R)')
xlabel('R [km]')
ylabel('SNR(R)') 

% Plot Noises

Rvector = [R_200;R_1;R_2;R_3d;R_3u;R_4];
Sigma_shs_vector = [Sigma_shs_200*B;Sigma_shs_1*B;Sigma_shs_2*B;Sigma_shs_3d*B;Sigma_shs_3u*B;Sigma_shs_4*B];
Sigma_shd_vector = [Sigma_shd*B;Sigma_shd*B;Sigma_shd*B;Sigma_shd*B;Sigma_shd*B;Sigma_shd*B];
Sigma_th_vector = [Sigma_th*B;Sigma_th*B;Sigma_th*B;Sigma_th*B;Sigma_th*B;Sigma_th*B];
figure
semilogy(Rvector,Sigma_shs_vector,':r*');
title('N(V^2)')
xlabel('R [km]')
ylabel('N(V^2)') 
hold on
semilogy(Rvector,Sigma_shd_vector,':b*');
hold on
semilogy(Rvector,Sigma_th_vector,':g*');
legend('shs','shd','th')

% 6. Assess the approximate laser-radar maximum range (SNR(Rmax)=1)

P_max = sqrt(Sigma_th*B)/(Rv*L);
R_20km = 0.2:0.1:20;

for i=1:199
   p_r(i) = power_return(R_20km(i), R_PBL, Alpha_aer, Alpha_mol, Beta_aer, Beta_mol, K)
   Sigma_shs_r(i) = 2*q*GT^2*F*M^2*Rio*(p_r(i)+Pback)*L; % [V^2/Hz]
   Sigma_v_r(i) = Sigma_shs_r(i)+Sigma_shd+Sigma_th;
   N_r(i) = Sigma_v_r(i)*B; % [V^2]
   SNR_r(i) = (Rv*L*p_r(i))/(sqrt(N_r(i)));
   SNR_r_db(i) = 20*log10(SNR_r(i));
end

in = find(SNR_r_db >= 1);
in_sol = max(in);
R_max = R_20km(in_sol)


%syms Rm K Beta_mol Alpha_aer Alpha_mol R_PBL P_max
%eqn = (K/Rm^2)*(Beta_mol)*exp(-2*(Alpha_aer+Alpha_mol)*R_PBL)*exp(-2*Alpha_mol*(Rm-R_PBL)) == P_max;
%Rmax = solve(eqn,Rm)

% 7. How many pulses are needed to integrate in order to ensure a SNRv (voltage signal-to-noise ratio) of 40 dB at 3-km range? What is the resulting observation time of the lidar instrument?

SNR_40dB = 10^(40/20)
Npulses = ceil(SNR_40dB^2/SNR_3d^2)
freq = 20; % [Hz]
Tobs = (1/20)*Npulses

% 8. Now, consider a Raman system of similar specs. If for Raman systems the return signal is typically 3 orders of magnitude lower than for their elastic system counterparts, discuss on the feasibility of day-time operation.

R_raman = linspace(0.2,4);
P_raman = 1e-3*(K./R_raman.^2).*(Beta_aer+Beta_mol).*exp(-2*(Alpha_aer+Alpha_mol).*R_raman); % [W]

figure
semilogy(R_raman,P_raman,'-');
title('Return Raman Power Praman(R)')
xlabel('R [km]')
ylabel('Praman(R) [W]') 
hold on
semilogy(Rvector,Pvector,':r*');
hold on
semilogy(R_raman,ones(size(R_raman))*3e-6,'g');
legend('Raman','Elastic','Day-time background')

% 9. Compute the photodiode NEP and its quantum efficiency (eta)

h = 6.6262e-34;
Lambda = 532e-9;

NEP = (sqrt(Sigma_shd/GT^2))/(Rio*M) % [W/sqrt(Hz)]
Eta = (Rio*h*c)/(q*Lambda)

% 10 Compute the system NEP (NEPs)

NEPs = (sqrt(Sigma_shd+Sigma_th))/(Rv_net) % [W/sqrt(Hz)]

