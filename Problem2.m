%Problem 2
close all;
clear all;
% 1. Determine the system constant, K [W·km3]

c = 2.99793e8; % [m/s]
E = 300e-3; % [J]
d_primary = 0.5; % [m]
R = (d_primary/2)^2-(0.06858/2)^2; % [m]
Ar = pi*R; % [m^2]
K_1 = ((E*c)/2)*Ar; % [Wm^3]
K = K_1*1e-9 % [Wkm^3]


% 2. Estimate the received background power for both the elastic and Raman channels (Pback,0 and Pback,R, respectively) under (see “day-time/night-time” parameter) operation.

Lsun = 3e-6; % [Wcm-2nm-1sr-1]
Ar_1 = Ar*1e4; % [cm^2]
rd = 3e-3/2; % [m]
f = 2; % [m]
FOV = atan(rd/f); % [rad]
Delta_Omega = pi*(sin(FOV))^2; % [sr]
Delta_Lambda_elastic = 10; %[nm]
Delta_Lambda_raman = 0.5; %[nm]
Pback_elastic = Lsun*Ar_1*Delta_Omega*Delta_Lambda_elastic % [W]
Pback_raman = Lsun*Ar_1*Delta_Omega*Delta_Lambda_raman % [W]

% 3. Plot both elastic- and Raman-return powers (P0(R) and PR(R), respectively). Superimpose Pback,0 and Pback,R plots from question 2 results.

R_15km = 0:0.1:15;
R_PBL = 3; % [km]

for i=1:151
   Alpha_mol = 1.2569e-2 - 7.7599e-4*R_15km(i); % [km-1]
   Alpha_aer = 1; % [km-1]
   Beta_aer = Alpha_aer/25; % [km-1]
   Beta_mol = Alpha_mol/(8*pi/3); % [km-1]
   p_r(i) = power_return(R_15km(i), R_PBL, Alpha_aer, Alpha_mol, Beta_aer, Beta_mol, K)
   NR = 2.1145e34 - 2.0022e33*R_15km(i) + 5.4585e31*R_15km(i)^2;
   RCross_section = 3.71e-41; % [km2sr-1]
   Lambda_r  = 607.4e-9; % [m]
   Lambda_0  = 532e-9; % [m]
   Alpha_mol_int = 1.2569e-2*R_15km(i) - 7.7599e-4*R_15km(i)*R_15km(i)/2;
   
   Alpha_raman_aer = ((Lambda_r/Lambda_0)^(-1.8))*Alpha_aer;
   Alpha_raman_mol = 7.3219e-3 - 4.5204e-4*R_15km(i); % [km-1]
   Alpha_raman_mol_int = 7.3219e-3*R_15km(i) - 4.5204e-4*R_15km(i)*R_15km(i)/2;
   p_r_raman(i) = (K/R_15km(i))*(NR*RCross_section)*exp(-2 * (Alpha_aer*R_15km(i) + Alpha_mol_int*+Alpha_raman_aer* R_15km(i)+Alpha_raman_mol_int) )
end



figure
semilogy(R_15km,p_r,'r');
hold on
semilogy(R_15km,p_r_raman,'b');
hold on
semilogy(R_15km,ones(size(R_15km))*Pback_elastic,'c');
hold on
semilogy(R_15km,ones(size(R_15km))*Pback_raman,'y');
title('Return Power P(R)')
xlabel('R [km]')
ylabel('P(R) [W]') 
legend('Elastic','Raman','Pback Elastic','Pback Raman')


% 4. Compute, for the elastic and Raman channels, receiver-chain voltage responsivities (Rv,0 and Rv,R, respectively), and net voltage responsivities (i.e., including spectral optical losses; Rv,0’ and Rv,R’, respectively).

% 5. a) Assuming analog detection, plot the elastic range-dependent signal-to-noise ratio, SNR0(R), at the output of the receiver chain (i.e., voltage ratio) and related shot photo-induced, shot-dark and thermal variances