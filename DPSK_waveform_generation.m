% DPSK waveform generation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RZ signal
% ===============================================================
% RZ_sig  = [];
% CLK_sig = [];

% sin wave

ii = 1:patt_len*bit_d;

CLK_sig = 0.5 * (sin(2 * pi * (ii' / bit_d) - pi / 2) + 1);

%{
%sech-shaped pulse
fwhm=bit_d/10;           %33[%]の形状
CLK = sech(2*log(1+sqrt(2))/fwhm*(-20:19));
CLK = CLK-CLK(:,1);
CLK_sig = CLK;

for ii=1:1:patt_len-1
   CLK_sig = horzcat(CLK_sig, CLK);
end
CLK_sig=CLK_sig';
%}
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

% input_pro1   = 10^(input_level_pro1 / 10) * 10^(-3) * ones(dat_len,1);
% input_pro2   = 10^(input_level_pro2 / 10) * 10^(-3) * ones(dat_len,1);

% NRZ signal
% ===============================================================

% Amp_sig    = ones(dat_len,1) * 10^(input_level_sig / 10) * 10^(-3);

%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NRZ_sig  = zeros(bit_d,1);

%t_rise = 0.8;  % rise time duration, ratio to bit-slot, should be less than 1.0 (1.0 => 0.999)
%t_fall = 0.8;  % fall time duration, ratio to bit-slot, should be less than 1.0 (1.0 => 0.999)

t_rise = 0.8;  % rise time duration, ratio to bit-slot, should be less than 1.0 (1.0 => 0.999)
t_fall = 0.8;  % fall time duration, ratio to bit-slot, should be less than 1.0 (1.0 => 0.999)

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
    %NRZ_sig?ｽﾌサ?ｽC?ｽY?ｽ?ｽbit_d?ｽ?ｽ?ｽg?ｽ?ｽ?ｽq?ｽA?ｽp?ｽ^?ｽ[?ｽ?ｽ?ｽ?ｽ1?ｽﾈゑｿｽtmp_NRZ_sig?ｽ?ｽﾇ会ｿｽ
    NRZ_sig  = [NRZ_sig;zeros(bit_d,1)] + pattern(1,i) * [zeros((i-1)*bit_d,1);tmp_NRZ_sig];

end
%bit_d?ｽﾌ費ｿｽ?ｽ?ｽ?ｽ?ｽ?ｽN?ｽ_?ｽ?ｽNRZ_sig?ｽz?ｽ?ｽﾌサ?ｽC?ｽY?ｽﾏ更
NRZ_sig = NRZ_sig(bit_d/2+1:bit_d/2+bit_d*patt_len);


dat_len = size(NRZ_sig,1);

NRZ_sig(dat_len,1) = NRZ_sig(dat_len-1,1);
NRZ_sig(dat_len/2,1) = NRZ_sig(dat_len/2-1,1);

E_input_sig   = sqrt(Amp_sig) .* exp(Im*pi*NRZ_sig);
P_input_sig   = Amp_sig;

Phase_input_sig = pi * NRZ_sig;


%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probe (Clock or CW)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%input_pu1   = 10^(input_level_pu1 / 10) * 10^(-3) * CLK_sig;
%input_pu2   = 10^(input_level_pu2 / 10) * 10^(-3) * CLK_sig;

input_pu1   = 10^(input_level_pu1 / 10) * 10^(-3) * ones(dat_len,1);
input_pu2   = 10^(input_level_pu2 / 10) * 10^(-3) * ones(dat_len,1);

t = 0:del_t:del_t*(dat_len-1);
% 
% 
% figure
% subplot(2,1,1)
% plot(t*10^9,P_input_sig*10^3);
% xlabel('Time[ns]');
% ylabel('Intensity[mW]')
% xlim([10*10^(-1)-10*10^(-1),15*10^(-1)-10*10^(-1)]);
% subplot(2,1,2)
% plot((t-10*10^(-10))*10^9,Phase_input_sig);
% xlabel('Time[ns]');
% ylabel('Phase[rad]')
% xlim([10*10^(-1)-10*10^(-1),15*10^(-1)-10*10^(-1)]);

