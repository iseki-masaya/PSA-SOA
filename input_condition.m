%%%%%%%%%%%%%%%%%%%%%%%%
%input condition
%%%%%%%%%%%%%%%%%%%%%%%%

I_current           = 500*10^(-3);			%injected current [A]

lam_sig             = 1.550e-6;         	%wavelength of signal [m]
lam_pu1             = 1.5495e-6;			%wavelength of pump light 1 [m]
lam_pu2             = 1.5505e-6;			%wavelength of pump light 2 [m]

del_lam             = (abs(lam_sig-lam_pu1)+abs(lam_pu2-lam_sig))/2;

eta_fwm             = 10^(-2.1164 * del_lam * 1e9 /10);    %from experiments
eta_hosei           = 1100*(del_lam*1e9)^2+1300*(del_lam*1e9)-700;

E_sig               = (h * c) / lam_sig;		%photon energy of signal light [J]
E_pu1               = (h * c) / lam_pu1;		%photon energy of pump light 1 [J]
E_pu2               = (h * c) / lam_pu2;		%photon energy of pump light 2 [J]

del_nyu             = 0.4928649e12;     %bandwidth of ASE [Hz]

pattern_gene_prbs5

patt_len            = size(pattern,2);			%pattern length
dat_len             = bit_d * patt_len;			%data length

%SOAÇÃì¸óÕëπ1.4dBÅCèoóÕëπ1.6dB
input_level_sig     = -3;			%input power of signal light [dBm]
input_level_pu1     = 15;			%input power of pump light 1 [dBm]
input_level_pu2     = 15;			%input power of pump light 2 [dBm]

beta                = 2*pi*c*(1/lam_pu1+1/lam_pu2-2/lam_sig)*del_t;	%phase mismatch

gamma               =(16000 + eta_hosei) * eta_fwm  * 2000e-6 / div_n; %é¿å±Ç…ÇÊÇÈîgí∑àÀë∂ÇéÊÇËì¸ÇÍÇΩåãâ 
disp('previous gamma')
disp(gamma)
% gamma               = 2.25e5 * 2000e-6 / div_n;
% disp('gamma')
% disp(gamma)