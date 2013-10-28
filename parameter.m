%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%?p?????[?^????
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h               = 6.6260693e-34;            %?v?????N????[J/s]
c               = 2.99792458e8;             %????[m/s]
Im              = sqrt(-1);                 %???f??
%gam             = 0.3;                      %???????????W??
gam             = 0.14;
C               = 2e-5;                     %???R???o???W??
%carrier_den_tr  = 0.7e24;
carrier_den_tr  = 1.1e24;                   %?????L?????A???x[m^(-3)]
%carrier_den_tr  = 0.9e24;                   %?????L?????A???x[m^(-3)]
%lam_tr          = 1.587415e-6;                 %?????g??[m]
lam_tr          = 1.585595e-6;                 %?????g??[m]
%lam_tr          = 1.605e-6;                 %?????g??[m]
%lam_tr          = 1.652e-6;                 %?????g??[m]
index_group     = 4;                        %?Q??????
index_ini       = 3.5;                      %??????????
% diff_ref_index    = -1.2e-26;		%differencial index (bulk)
diff_ref_index      = -1.8e-27;		%differencial index (MQW)
group_velocity  = c / index_group;          %?Q???x[m/s]
q_unit          = 1.60217646e-19;           %?d???f??

% recom_rate1     = 4e8;                    %?????????[?g[s^(-1)]
% recom_rate2     = 4e-16;                  %?????????[?g[m^3/s]
% recom_rate3     = 1e-40;                  %?????????[?g[m^6/s]

%recom_rate1     = 1.0e8;                    %?????????[?g[s^(-1)]
%recom_rate2     = 2.5e-17;                  %?????????[?g[m^3/s]
%recom_rate3     = 9.4e-41;                  %?????????[?g[m^6/s]

recom_rate1         = 3*1.0e8;		%recombination rate 1[s^(-1)]
recom_rate2         = 3*2.5e-17;		%recombination rate 2[m^3/s]
recom_rate3         = 5*9.4e-41;		%recombination rate 3[m^6/s]

% recom_rate1     = 13.5e8;                    %?????????[?g[s^(-1)]
% recom_rate2     = 56e-17;                  %?????????[?g[m^3/s]
% recom_rate3     = 15e-41;                  %?????????[?g[m^6/s]


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

%gain_compress   = 0;                  %?????`???????k???q[m^3]
gain_compress   = 1.3e-23;                  %?????`???????k???q[m^3](?????l)
%gain_compress   = 1.3e-23;                  %?????`???????k???q[m^3]

%loss_in_active  = 140e2;                    %?????w????????[m^(-1)]
%loss_in_active  = 2.0e3;                    %?????w????????[m^(-1)]

k_0                 = 9000;         %active layer loss constant 0 [m^-1]
k_1                 = 2250e-24;     %active layer loss constant 1 [m^-1]

loss_in_clad    = 20e2;                    %?N???b?h?w????????[m^(-1)]
%loss_in_clad    = 2.0e2;                    %?N???b?h?w????????[m^(-1)]
scatt_loss      = 1.0e2;                    %?U??????[m^(-1)]

Ref             = 0.00001;                   %?[????????????

quantum_eff     = 0.5;                      %???q????

Be              = 10e9;                     %?d?C?G????????

%carrier_den_0   = 5.669e24;                 %?????L?????A???x[m^(-3)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%?M?G??
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
boltzmann       = 1.3806523e-23;            %?{???c?}??????[J/K]
term            = 300;                      %???x[K]
Res             = 50;                       %?????G?????R[Ohm]
%?M?G??
N_th            = (4 * boltzmann * term * Be) / Res;
sigma_th        = sqrt(N_th);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%?f?o?C?X?p?????[?^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% length          = 900e-6;                   %?f?o?C?X??[m]
% 
% volume          = length * A_cross;         %?????w????[m^3]

% del_z           = length / div_n;           %??????????????[m]
% del_t           = del_z / group_velocity;   %????????[s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%?????M???p?????[?^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I_current       = 400*10^(-3);              %?????d??[A]
%I_current       = 300*10^(-3);              %?????d??[A]
%lam_sig         = 1.5557e-6;                 %??????1?g??[m]
%lam_pro1        = 1.5541e-6;                 %?v???[?u???g??[m]
%lam_pro2        = 1.5541e-6;                 %?v???[?u???g??[m]
%lam_sig         = 1.555e-6;                 %??????1?g??[m]
%lam_sig2         = 1.555e-6;                 %??????2?g??[m]
%lam_pro1        = 1.550e-6;                 %?v???[?u???g??[m]
%lam_pro2        = 1.550e-6;                 %?v???[?u???g??[m]
%E_sig           = (h * c) / lam_sig;       %??????1???G?l???M?[
%E_sig2           = (h * c) / lam_sig2;       %??????2???G?l???M?[
%E_pro1          = (h * c) / lam_pro1;        %?v???[?u?????G?l???M?[
%E_pro2          = (h * c) / lam_pro2;        %?v???[?u?????G?l???M?[

%------------ device parameter------------

wide                = 1.3e-6;		%width of active layer [m]
depth               = 0.1e-6;		%depth of active layer [m]
A_cross             = wide * depth * 0.7;  	%cross section of active layer [m^2](well:varrier = 7:3)
length              = 2000e-6;		%device length [m]
volume              = length * A_cross;	%volume of active layer [m^3]

speed           = 40e9;                     %?????M?????r?b?g???[?g
bit_d           = 2;                       %1?r?b?g?f?[?^???????? 10
del_t           = 1 / (speed * bit_d);      %????????[s]
del_z           = del_t * group_velocity;   %??????????????[m]
div_n               = floor(length / del_z);	%devided number of SOA

fs              = 1/del_t;