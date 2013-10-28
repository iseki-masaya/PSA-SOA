%%%%%%%%%%%%%%%%%%%%%%%%
%input condition
%%%%%%%%%%%%%%%%%%%%%%%%

I_current           = 500*10^(-3);			%injected current [A]

lam_sig             = 1.550e-6;         	%wavelength of signal [m]
lam_pu1             = 1.5495e-6;			%wavelength of pump light 1 [m]
lam_pu2             = 1.5505e-6;			%wavelength of pump light 2 [m]

%{
wavelength spacing = 1nm
lam_sig             = 1.5505e-6;			%wavelength of signal [m]
lam_pu1             = 1.5500e-6;			%wavelength of pump light 1 [m]
lam_pu2             = 1.5510e-6;			%wavelength of pump light 2 [m]
%}
%del_fre             = (abs(1/lam_sig-1/lam_pu2)+abs(1/lam_sig-1/lam_pu1))*c/2*1e-12;
%eta_fwm             = -8.778 * log(del_lam) -1.8992;
%eta_fwm             = -14.998*del_fre+0.5928;
%eta_fwm             = 10 ^ (eta_fwm/10);

del_lam             = (abs(lam_sig-lam_pu1)+abs(lam_pu2-lam_sig))/2;

eta_fwm             = 10^(-2.1164 * del_lam * 1e9 /10);    %from experiments
eta_hosei           = 1100*(del_lam*1e9)^2+1300*(del_lam*1e9)-700;

%{
lam_sig             = 1.550e-6;			%wavelength of signal [m]
lam_pu1             = 1.540e-6;			%wavelength of pump light 1 [m]
lam_pu2             = 1.560e-6;			%wavelength of pump light 2 [m]


lam_sig             = 1.550e-6;			%wavelength of signal [m]
lam_pu1             = 1.5475e-6;			%wavelength of pump light 1 [m]
lam_pu2             = 1.5525e-6;			%wavelength of pump light 2 [m]
%}
E_sig               = (h * c) / lam_sig;		%photon energy of signal light [J]
E_pu1               = (h * c) / lam_pu1;		%photon energy of pump light 1 [J]
E_pu2               = (h * c) / lam_pu2;		%photon energy of pump light 2 [J]

%speed               = 40e9;				%bit rate [bps]
%bit_d               = 40;				%sample number in one bit time slot
%del_t               = 1 / (speed * bit_d);		%time sampling interval [s]
%del_z               = del_t * group_velocity;	%space sampling interval [m]
%div_n               = floor(length / del_z);	%devided number of SOA

del_nyu             = 0.4928649e12;     %bandwidth of ASE [Hz]

%pattern_gene_impulse40G				%pattern generation 
%pattern_gene_step
pattern_gene_prbs15
%pattern_gene_prbs7
%pattern_gene_prbs6
%pattern_gene_prbs5
%pattern_gene_cw
%pattern_gene_psk
%pattern_gene_prbs10

patt_len            = size(pattern,2);			%pattern length
dat_len             = bit_d * patt_len;			%data length

%SOAの入力損1.4dB，出力損1.6dB
input_level_sig     = -3;			%input power of signal light [dBm]
input_level_pu1     = 15;			%input power of pump light 1 [dBm]
input_level_pu2     = 15;			%input power of pump light 2 [dBm]

beta                = 2*pi*c*(1/lam_pu1+1/lam_pu2-2/lam_sig)*del_t;	%phase mismatch

%gamma              = 155150 * 2000e-6 / div_n;	%nonlinear coefficient per devided unit [1/W]
%gamma               = 17136 * 2000e-6 / div_n;	%nonlinear coefficient per devided unit [1/W]
%gamma               = 5000 * 2000e-6 / div_n;	%nonlinear coefficient per devided unit [1/W]
%gamma              =  0 * 2000e-6 / div_n;
%gamma               =9559 * 2000e-6 / div_n * eta_fwm;    %参考文献のQD-SOAでかつ波長依存性を考慮するとこの値がぴったり
gamma               =(16000 + eta_hosei) * eta_fwm  * 2000e-6 / div_n; %実験による波長依存を取り入れた結果
