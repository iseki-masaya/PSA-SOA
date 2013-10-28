% waveform generation for dpsk
% Amplitude
h                   = 6.6260693e-34;	%planc constant [J/s]
c                   = 2.99792458e8;		%velocity of light [m/s]
%pi                 = 3.141592653589793	%pi
Im                  = sqrt(-1);		%imaginary number
gam                 = 0.14;			%confinement factor
C                   = 2e-5;			%spontaneous emission coefficient
carrier_den_tr      = 1.1e24;		%carrier density at transparency [m^(-3)]
% lam_tr            = 1.587415e-6;	%wavelength at transparency [m]
index_group         = 4;			%group index
index_ini           = 3.5;			%initial index
% diff_ref_index    = -1.2e-26;		%differencial index (bulk)
diff_ref_index      = -1.8e-27;		%differencial index (MQW)
group_velocity      = c / index_group;	%group velocity [m/s]
q_unit              = 1.60217646e-19;	%electron charge [C]

recom_rate1         = 3*1.0e8;		%recombination rate 1[s^(-1)]
recom_rate2         = 3*2.5e-17;		%recombination rate 2[m^3/s]
recom_rate3         = 5*9.4e-41;		%recombination rate 3[m^6/s]

%gain_const0        = 2.5e-20;		%gain constant 0 [m^2]
%gain_const1        = 0.074e20;		%gain constant 1[m^(-3)]
%gain_const2        = 3.0e-32;		%gain constant 2[m^4]
%gain_const3        = 3.155e25;		%gain constant 3[m^(-4)]

gain_const3         = 2.07e25;		%gain constant 3[m^(-4)]

%------------------ peak-gain -----------------
a_g1                = 2.6565436412E-20;                   
b_g1                = 1.3216883297E+03;

a_g2                = 5.1998178058E-21;
b_g2                = 4.5215586037E+04;

cg                  = 3.66e6;

%----------- polynomial coefficient -----------
a_p1                = 2.2183493359E-06;
b_p1                = -1.3202457757E+16;

a_p2                = 3.0724794162E-07;
b_p2                = 4.1892792542E+18;

cp                  = 2.88e34;

%------------ peak-gain wavelength------------
a_l1                = -6.4831485654E-32;
b_l1                = 1.6878978292E-06;

a_l2                = -7.1874120023E-33;
b_l2                = 1.5743719969E-06;

cl                  = 1.78e-17;

gain_compress       = 1.3e-23;		%nonlinear gain compression constant[m^3]

loss_in_active      = 140e2;		%loss coefficient in active layer [m^(-1)]
loss_in_clad        = 20e2;			%loss coefficient in clad layer [m^(-1)]
scatt_loss          = 1.0e2;		%scattering loss [m^(-1)]

Ref                 = 0.00001;		%edge reflectivity

quantum_eff         = 0.5;			%quantum efficiency

%------------ device parameter------------

wide                = 1.3e-6;		%width of active layer [m]
depth               = 0.1e-6;		%depth of active layer [m]
A_cross             = wide * depth * 0.7;  	%cross section of active layer [m^2](well:varrier = 7:3)
length              = 2000e-6;		%device length [m]
volume              = length * A_cross;	%volume of active layer [m^3]

%------------ FWM parameter------------

%beta               = 0;			%phase mismatching

I_current           = 300*10^(-3);			%injected current [A]

lam_sig             = 1.55025e-6;			%wavelength of signal [m]
lam_pu1             = 1.55000e-6;			%wavelength of pump light 1 [m]
lam_pu2             = 1.55050e-6;			%wavelength of pump light 2 [m]

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

speed               = 40e9;				%bit rate [bps]
bit_d               = 40;				%sample number in one bit time slot
del_t               = 1 / (speed * bit_d);		%time sampling interval [s]
del_z               = del_t * group_velocity;	%space sampling interval [m]
div_n               = floor(length / del_z);	%devided number of SOA

pattern_gene_impulse40G				%pattern generation 
%pattern_gene_step
%pattern_gene_prbs6
%pattern_gene_prbs4
%pattern_gene_cw

patt_len            = size(pattern,2);			%pattern length
dat_len             = bit_d * patt_len;			%data length

input_level_sig     = -3;			%input power of signal light [dBm]
input_level_pu1     = 10;			%input power of pump light 1 [dBm]
input_level_pu2     = 10;			%input power of pump light 2 [dBm]

beta                = 2*pi*c*(1/lam_pu1+1/lam_pu2-2/lam_sig)*del_t;	%phase mismatch

%gamma              = 155150 * 2000e-6 / div_n;	%nonlinear coefficient per devided unit [1/W]
%gamma               = 17136 * 2000e-6 / div_n;	%nonlinear coefficient per devided unit [1/W]
gamma               = 11500 * 2000e-6 / div_n;	%nonlinear coefficient per devided unit [1/W]
%gamma              =  0 * 2000e-6 / div_n;

ii = 1:patt_len*bit_d;

CLK_sig = 0.5 * (sin(2 * pi * (ii' / bit_d) - pi / 2) + 1);
RZ_sig = CLK_sig;

% for i = 1:1:patt_len
% 
%     tmp_RZ_sig  = [];
%     tmp_CLK_sig = [];
% 
%     for ii = 1:1:bit_d
%         tmp_sig_pow = 0.5 * (sin(2 * pi * (ii / bit_d) - pi / 2) + 1);
%         tmp_CLK_pow = 0.5 * (sin(2 * pi * (ii / bit_d) - pi / 2) + 1);
%         tmp_RZ_sig  =  [tmp_RZ_sig;tmp_sig_pow];
%         tmp_CLK_sig =  [tmp_CLK_sig;tmp_CLK_pow];
%     end
% 
%     RZ_sig  = [RZ_sig;tmp_RZ_sig];
%     CLK_sig = [CLK_sig;tmp_CLK_sig];
% 
% end

% RZ_sig = [zeros(bit_d/2,1);RZ_sig;zeros(bit_d/2,1)];

clear tmp_RZ_sig tmp_sig_pow tmp_CLK_sig tmp_CLK_pow i ii;

CLK_delay_t = 6.25e-12;
CLK_delay_n = floor(CLK_delay_t/del_t);

tmp_CLK_sig = CLK_sig(1:dat_len-bit_d,1);

CLK_sig = [zeros(CLK_delay_n,1);tmp_CLK_sig;zeros(bit_d-CLK_delay_n,1)];

Amp_sig   = RZ_sig * 10^(input_level_sig / 10) * 10^(-3);

% Phase

NRZ_sig  = zeros(bit_d,1);

t_rise = 0;  % rise time duration MAX=0.999
t_fall = 0;  % fall time duration MAX=0.999

bit_d_rise = floor(bit_d*t_rise*0.5);
bit_d_fall = floor(bit_d*t_fall*0.5);

tmp_NRZ_sig  = zeros(2*bit_d,1);
tmp_NRZ_sig(bit_d/2+1:bit_d*3/2,1) = 1;
tmp_NRZ_sig(bit_d/2,1) = 0.5;
tmp_NRZ_sig(bit_d*3/2,1) = 0.5;

% ii = 1:1:bit_d_rise;
%     tmp_NRZ_sig(bit_d/2 + ii, 1)  = 0.5 * (sin(2 * pi * (ii / (2*t_rise*bit_d))) + 1);
%     tmp_NRZ_sig(bit_d/2 - ii, 1)  = 0.5 * (sin(2 * pi * (-ii / (2*t_rise*bit_d))) + 1);

for ii = 1:1:bit_d_rise
    tmp_NRZ_sig(bit_d/2 + ii, 1)  = 0.5 * (sin(2 * pi * (ii / (2*t_rise*bit_d))) + 1);
    tmp_NRZ_sig(bit_d/2 - ii, 1)  = 0.5 * (sin(2 * pi * (-ii / (2*t_rise*bit_d))) + 1);
end

for ii = 1:1:bit_d_fall
    tmp_NRZ_sig(bit_d*3/2 + ii, 1)  = 0.5 * (sin(2 * pi * (-ii / (2*t_fall*bit_d))) + 1);
    tmp_NRZ_sig(bit_d*3/2 - ii, 1)  = 0.5 * (sin(2 * pi * (ii / (2*t_fall*bit_d))) + 1);
end

for i = 1:1:patt_len
    %NRZ_sig�̃T�C�Y��bit_d���g���q�A�p�^�[����1�Ȃ�tmp_NRZ_sig��ǉ�
    NRZ_sig  = [NRZ_sig;zeros(bit_d,1)] + pattern(1,i) * [zeros((i-1)*bit_d,1);tmp_NRZ_sig];

end
%bit_d�̔������N�_��NRZ_sig�z��̃T�C�Y�ύX
NRZ_sig = NRZ_sig(bit_d/2+1:bit_d/2+bit_d*patt_len);

dat_len = size(NRZ_sig,1);
E_input_sig   = sqrt(Amp_sig) .* exp(Im*pi*NRZ_sig);
P_input_sig   = Amp_sig;
Phase_input_sig = pi * NRZ_sig;

% Probe (Clock or CW)

% input_pro1   = 10^(input_level_pro1 / 10) * 10^(-3) * CLK_sig;
% input_pro2   = 10^(input_level_pro2 / 10) * 10^(-3) * CLK_sig;

input_pu1   = 10^(input_level_pu1 / 10) * 10^(-3) * ones(dat_len,1);
input_pu2   = 10^(input_level_pu2 / 10) * 10^(-3) * ones(dat_len,1);

t = 0:del_t:del_t*(dat_len-1);
%plot(t,angle(E_input_sig));
plot(t,E_input_sig);