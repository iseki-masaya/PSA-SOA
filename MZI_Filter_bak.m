function [E_out, Power, Phase,Q_map,sigma_map] = MZI_Filter(Q_map,...
    sigma_map,i_alpha,i_i_s,input_waveform,dt, input_phase, D, S, L, wavelength, center_wavelength,if_Gaussnoise,if_sa,if_fluc)
%function [complex_waveform, Power, Phase] = MZI_Filter(input_waveform,dt, input_phase, D, S, L, wavelength, center_wavelength)
close all
parameter;
% input_condition;

if_fig = 0;
% data size verify
% ------------------------------------------

dat_size_P_r = size(input_waveform,1);
dat_size_P_c = size(input_waveform,2);

if dat_size_P_r == 1
    if dat_size_P_c == 1
        error('signal data is invalid !!');
    end
    dat_size_P = dat_size_P_c;
    input_waveform = input_waveform';
elseif dat_size_P_c == 1
    dat_size_P = dat_size_P_r;
end

dat_size_phase_r = size(input_phase,1);
dat_size_phase_c = size(input_phase,2);

if dat_size_phase_r == 1
    if dat_size_phase_c == 1
        error('signal data is invalid !!');
    end
    dat_size_phase = dat_size_phase_c;
    input_phase = input_phase';
elseif dat_size_phase_c == 1
    dat_size_phase = dat_size_phase_r;
end

if dat_size_P ~= dat_size_phase
    error('signal data is invalid !!');
end

%??�ǉ�
%t = 0:del_t:del_t*(dat_size_P-1);

t = get_time(dat_size_P)';
n = 1:1:dat_size_P;

% -------------------------------------------

E_in = sqrt(input_waveform) .* exp(1i*input_phase);
% E_in = signal_gain(E_in,P_i);
% [E_in, ~, ~] = SPM(input_waveform,...
%     del_t,input_phase,D_SSMF,S_SSMF,100e3,1550e-9, 1550e-9); 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%??�??�??�x??�G??�??�%%%%%
if(if_Gaussnoise)
noise_size = size(t,1);

[N_ssp,N_spsp,N_RIN,~,N_spshot,N_shot_p] = gaussian_noise();
noise_ssp = (randn(noise_size,1) .* sqrt(N_ssp));    %??�d??�??�??�U??�??�
noise_spsp = (randn(noise_size,1) .* sqrt(N_spsp));      %??�d??�??�??�U??�??�
noise_RIN =  (randn(noise_size,1) .* sqrt(N_RIN)); 
noise_spshot =  (randn(noise_size,1) .* sqrt(N_spshot)); 
noise_shot_p =  (randn(noise_size,1) .* sqrt(N_shot_p));

noise_phase1 = exp(1i .* rand(noise_size,1) * pi());
noise_phase2 = exp(1i .* rand(noise_size,1) * pi());
noise_phase3 = exp(1i .* rand(noise_size,1) * pi());
noise_phase4 = exp(1i .* rand(noise_size,1) * pi());
% noise_phase1 = 1;
% noise_phase2 = 1;
% noise_phase3 = 1;
% noise_phase4 = 1;
SNR = 10*log10(mean(abs(E_in).^2)/(N_ssp+N_spsp+N_RIN+N_spshot+N_shot_p));
E_in = awgn(E_in,21.6,'measured');

noise_t = (noise_ssp.*noise_phase1 + noise_spsp.*noise_phase2+noise_RIN+...
    noise_spshot.*noise_phase3+noise_shot_p.*noise_phase4);
% figure(3)
% plot(t*10^9,abs(noise_t).^2*10^3)
% xlabel('time[ns]');
% ylabel('Intensity[mW]');
% noise_temp = [];
% 
% for i = 1:(noise_size + 1)/noise_div
%     for ii = 1:noise_div
%         noise_temp = [noise_temp;noise_t(i,1)];
%     end
% end
% 
% size(E_in)
% size(noise_temp )
% noise_t =noise_temp;
clear noise_temp

% E_in = E_in + noise_t;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%E_in =awgn(E_in_power,SNR,'measured').*exp(1i * E_in_phase);
%E_in =awgn(E_in,SNR,'measured');
%??�ʑ�??�G??�??�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(if_fluc)
noise_size = size(t,1);
phase_delta = pi()/4;
% for i_fluc = 160:1:1000
%     
% E_in = (E_in.*exp(-1i* phase_delta .*randn(noise_size,1)));
E_in = (E_in.*exp(-1i* phase_delta .*sin(2*pi()*n/500))');
% Generate_Constellation(E_in2,bit_d);
% i_fluc
% pause(0.5)
% if(mod(i_fluc,30) == 0)
% close all
% end
% end
%5000??�??�??�炢
end

tn = t * 10^9;
figure(3)
subplot(2,1,1)
plot(tn-5.4,abs(E_in).^2 * 10^3);
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Intensity[mW]');
subplot(2,1,2)
plot(tn-5.4,angledef(E_in));
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Phase[rad]');

% i_f=0.9;
% i_b=0.945;
%                
%         E_in = SESAM(E_in,i_f,i_b,1.2,0.6*1e-4);
% %         E_in = signal_gain(E_in,-10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% figure(4)
% subplot(2,1,1)
% plot(tn-5.4,abs(E_in).^2 * 10^3);
% xlim([0,0.5]);
% xlabel('time[ns]');
% ylabel('Intensity[mW]');
% subplot(2,1,2)
% plot(tn-5.4,angledef(E_in));
% xlim([0,0.5]);
% xlabel('time[ns]');
% ylabel('Phase[rad]');

sigma_phai_in = var(angle(E_in));
%F_in = fftshift(fft(E_in));
% F_in = fft(E_in);

fs = 1/dt;

df = fs/dat_size_P;

%n = 0:1:(dat_size_P-1);
n = 0:1:(size(t,2)-1);

%f = (n*df - fs/2)';
f = (n*df)';

f0 = 2.99792458e8/center_wavelength;

fc = 2.99792458e8/wavelength;
 
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MZI_matrix = 1/sqrt(2)*[1,-1i;-1i,1]
%MZI_matrix =[1,-1i;1i,1];
%
%E_bar_1 = E_in;
%E_cross_1 = 0;
%E_bar_2 = [];
%E_cross_2 = [];
%E_bar_3 = [];
%E_cross_3 = [];
%
%E_bar_2 = MZI_matrix(1,1) * E_bar_1 + pi_rev(E_cross_1);
%E_cross_2 =pi_rev(E_bar_1); + MZI_matrix(2,2) * E_cross_1;
%
%%1bit??�x??�??�??�??�??�??�
%E_bar_2 = circshift(E_bar_2, [bit_d, 0]);
%E_bar_2(1:bit_d,1)=zeros(bit_d,1);
%
%E_bar_3 = MZI_matrix(1,1) * E_bar_2 +pi_rev(E_cross_2);
%E_cross_3 = pi_rev(E_bar_2) + MZI_matrix(2,2) * E_cross_2;
%subplot(2,1,1); plot(t,E_in)
%subplot(2,1,2); plot(t,E_cross_3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MZI_matrix = 1/sqrt(2)*[1,-1i;-1i,1]
MZI_matrix =[1,-1i;-1i,1]/sqrt(2);
MZI_waveform = {E_in;E_in * 0};

% %??�f??�[??�^??�擾??�p
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(t,E_in)
% xlabel('t (second)')
% ylabel('Intensity')
% saveas(gcf,'E_in.jpg')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %1bit??�x??�??�??�??�??�??�
% E_delay = circshift(E_in, [bit_d, 0]);
% 
% %??�ʑ�??�Ή�]
% E_pi = phase_rev(E_in);
% 
% E_bar = (E_delay + E_pi)/2;
% E_cross = (-E_delay-E_in)/2;

%3dB??�J??�v??�??�1??�??�??�
MZI_bar = (MZI_matrix(1,1) * MZI_waveform{1,1} + MZI_matrix(1,2) * MZI_waveform{2,1});
MZI_cross = (MZI_matrix(2,1) * MZI_waveform{1,1} + MZI_matrix(2,2) * MZI_waveform{2,1});

MZI_waveform{1,1} = MZI_bar;
MZI_waveform{2,1} = MZI_cross;


%%%%%%%%%%%%%%%%%%%??�P??�ʕϊ�%%%%%%%%%%%%%%%%%%%%%%%%
tn = t * 10^9;
if(if_fig)
figure(1)
subplot(8,2,1);
plot(tn-5.4,abs(MZI_bar).^2 * 10^3);
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Intensity[mW]');
title('(a)')
% saveas(gcf,'1.jpg');


subplot(8,2,2);
plot(tn-5.4,abs(MZI_cross).^2 * 10^3);
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Intensity[mW]')
title('(i)')
% saveas(gcf,'2.jpg');

subplot(8,2,3);
plot(tn-5.4,angledef(MZI_bar));
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Phase[rad]');
title('(b)')
% saveas(gcf,'3.jpg');

% ylim([-4,4]);
subplot(8,2,4);
plot(tn-5.4,angledef(MZI_cross));
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Phase[rad]');
title('(j)')
% saveas(gcf,'4.jpg');
end
%bar??�??�1bit+??�Βx??�??�
MZI_bar = [zeros(bit_d,1);MZI_bar];
MZI_cross = [MZI_cross;zeros(bit_d,1)];


t = get_time(size(MZI_bar,1))';
tn = t * 10^9;
MZI_waveform{1,1} = MZI_bar;
MZI_waveform{2,1} = MZI_cross;

if(if_fig)
subplot(8,2,5);
plot(tn-5.4,abs(MZI_bar).^2 * 10^3);
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Intensity[mW]');
title('(c)')
% saveas(gcf,'5.jpg');

subplot(8,2,6);
plot(tn-5.4,abs(MZI_cross).^2 * 10^3);
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Intensity[mW]');
title('(k)')
% saveas(gcf,'6.jpg');

% ylim([-4,4]);
subplot(8,2,7);
plot(tn-5.4,angledef(MZI_bar));
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Phase[rad]');
title('(d)')
% saveas(gcf,'7.jpg');

% ylim([-4,4]);
subplot(8,2,8);
plot(tn-5.4,angledef(MZI_cross));
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Phase[rad]');
title('(l)')
end
% saveas(gcf,'8.jpg');

%3dB??�J??�v??�??�2??�??�??�
MZI_bar = (MZI_matrix(1,1) * MZI_waveform{1,1} + MZI_matrix(1,2) * MZI_waveform{2,1});
MZI_cross = (MZI_matrix(2,1) * MZI_waveform{1,1} + MZI_matrix(2,2) * MZI_waveform{2,1});

MZI_waveform{1,1} = MZI_bar;
MZI_waveform{2,1} = MZI_cross;

% figure(9)
% plot(tn,abs(MZI_bar).^2)
% xlabel('time[ns]');
% ylabel('Intensity[mW]');

%??�??�??�`??�}??�??�??�ɂ�??�?�?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%180cm^-1,100??�??�
if(if_sa)
% MZI_bar = SESAM(MZI_bar,0.7,0.945);
% MZI_cross = SESAM(MZI_cross,0.7,0.945);

% figure(10)
% plot(tn,abs(MZI_bar).^2 * 10^3)
% xlabel('time[ns]');
% ylabel('Intensity[mW]');
% for i = 0.7:0.01:1
%     i_f = i;
%     for ii = 0.7:0.01:1
%         i_b = ii;
% i_alpha = alpha_0(d);
% i_i_s = I_sat(d);
% 
% i_f=0.7;
% i_b=0.945;
% 
%         MZI_bar = SESAM(MZI_bar,i_f,i_b,0.04,2*10^(-6));
%         MZI_cross = SESAM(MZI_cross,i_f,i_b,0.04,2*10^(-6));
%         
%         MZI_bar = signal_gain(MZI_bar,-10);
%         MZI_cross = signal_gain(MZI_cross,-10);
% i_f=0.9;
% i_b=0.945;
%                
%         MZI_bar = SESAM(MZI_bar,i_f,i_b,1.2,0.6*1e-4);
%         MZI_cross = SESAM(MZI_cross,i_f,i_b,1.2,0.6*1e-4); 
% % 
i_f=0.9;
i_b=0.945;
               
        MZI_bar = SESAM(MZI_bar,i_f,i_b,1.2,0.6*1e-4);
        MZI_cross = SESAM(MZI_cross,i_f,i_b,1.2,0.6*1e-4); 
        MZI_bar = signal_gain(MZI_bar,-10);
        MZI_cross = signal_gain(MZI_cross,-10);

 i_f=0.7;
i_b=0.945;

        MZI_bar = SESAM(MZI_bar,i_f,i_b,0.04,2*10^(-6));
        MZI_cross = SESAM(MZI_cross,i_f,i_b,0.04,2*10^(-6));
        


% 
% figure(10)
% % plot(tran_nl)
% subplot(2,1,1)
% plot(temp_intensity_b*10^3,tran_nl);
% xlabel('Intensity[mW]');
% ylabel('transmittance');
% 
% subplot(2,1,2)
% plot(temp_intensity_b*10^3);
% xlim([10,15*160]);
% xlabel('index');
% ylabel('Intensity[mW]');
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MZI_bar = MZI_bar * sqrt(10^3);
% MZI_cross = MZI_cross * sqrt(10^3);  

I_comp = 1 * 10^(-8);
% I_s_max = max(abs(MZI_cross).^2)
I_s_max = 1000;


for i = 1:1:size(t,1)
    if abs(MZI_cross(i,1))^2 < I_comp;
        MZI_cross(i,1) = 0;
    end
    if abs(MZI_cross(i,1))^2 > I_s_max;
        MZI_cross(i,1) = I_s_max;
    end
end
    
for i = 1:1:size(t,1)
    if abs(MZI_bar(i,1))^2 < I_comp
        MZI_bar(i,1) = 0;
    end
    if abs(MZI_bar(i,1))^2 > I_s_max
        MZI_bar(i,1) = I_s_max;
    end
end
if(if_fig)
    
figure(1)
subplot(8,2,9);
plot(tn-5.4,abs(MZI_bar).^2 * 10^3);
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Intensity[mW]');
title('(e)')
% saveas(gcf,'9.jpg');

subplot(8,2,10);
plot(tn-5.4,abs(MZI_cross).^2 * 10^3);
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Intensity[mW]');
title('(m)')
% saveas(gcf,'10.jpg');

% ylim([-4,4]);
subplot(8,2,11);
plot(tn-5.4,angledef(MZI_bar));
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Phase[rad]');
title('(f)')
% saveas(gcf,'11.jpg');

% ylim([-4,4]);
subplot(8,2,12);
plot(tn-5.4,angledef(MZI_cross));
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Phase[rad]');
title('(n)')
% saveas(gcf,'12.jpg');
end

MZI_waveform{1,1} = MZI_bar;
MZI_waveform{2,1} = MZI_cross;

MZI_bar = (MZI_matrix(1,1) * MZI_waveform{1,1} + MZI_matrix(1,2) * MZI_waveform{2,1});
MZI_cross = (MZI_matrix(2,1) * MZI_waveform{1,1} + MZI_matrix(2,2) * MZI_waveform{2,1});


% MZI_matrix_inv =[1,1i;1i,1]/(sqrt(2));
% 
% %??�܂�Ԃ�1??�??�??�
% 
% MZI_waveform{1,1} = MZI_bar;
% MZI_waveform{2,1} = MZI_cross;
% 
% MZI_bar = (MZI_matrix_inv(1,1) * MZI_waveform{1,1} + MZI_matrix_inv(1,2) * MZI_waveform{2,1});
% MZI_cross = (MZI_matrix_inv(2,1) * MZI_waveform{1,1} + MZI_matrix_inv(2,2) * MZI_waveform{2,1});
% 
% %cross??�??�1bit??�x??�??�
% MZI_bar = [MZI_bar;zeros(bit_d,1)];
% MZI_cross = [zeros(bit_d,1);MZI_cross];
% 
% %??�܂�Ԃ�2??�??�??�
% MZI_waveform{1,1} = MZI_bar;
% MZI_waveform{2,1} = MZI_cross;
% 
% MZI_bar = (MZI_matrix_inv(1,1) * MZI_waveform{1,1} + MZI_matrix_inv(1,2) * MZI_waveform{2,1});
% MZI_cross = (MZI_matrix_inv(2,1) * MZI_waveform{1,1} + MZI_matrix_inv(2,2) * MZI_waveform{2,1});

MZI_bar = MZI_bar * exp(1i * pi());

%??�??�0??�̃f??�[??�^??�̗}??�??�??�??�??�??�
for i = 1:1:size(t,1)
    if abs(MZI_cross(i,1))^2 < I_comp;
        MZI_cross(i,1) = 0;
    end
end
    
for i = 1:1:size(t,1)
    if abs(MZI_bar(i,1))^2 < I_comp
        MZI_bar(i,1) = 0;
    end
end



%%%%%%%%%%%??�ŏ�??�ƍŌ�??�0??�??�??�?�?%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:1:bit_d
%     MZI_bar(1,:) = [];
    MZI_bar(end,:) = [];
    MZI_cross(1,:) = [];
%     MZI_cross(end,:) = [];
end

temp_size = size(MZI_bar,1);
t = get_time(temp_size) ;
tn = t * 10^9;

% MZI_bar = MZI_bar * exp(1i * pi());
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = get_time(size(MZI_bar,1))';
tn = t * 10^9;
if(if_fig)
subplot(8,2,13);
plot(tn-5.4,abs(MZI_bar).^2 * 10^3);
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Intensity[mW]');
title('(g)')
% saveas(gcf,'13.jpg');

subplot(8,2,14);
plot(tn-5.4,abs(MZI_cross).^2 * 10^3);
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Intensity[mW]');
title('(o)')
% saveas(gcf,'14.jpg');

% ylim([-4,4]);
subplot(8,2,15);
%plot(t*10^9,unwrap(angle(MZI_bar)));
plot(tn-5.4,angle(MZI_bar));
ylim([-4,4]);
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Phase[rad]');
title('(h)')
% saveas(gcf,'15.jpg');

% ylim([-4,4]);
subplot(8,2,16);
plot(tn-5.4,angle(MZI_cross));
ylim([-4,4]);
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Phase[rad]');
title('(p)')
% saveas(gcf,'16.jpg');
end

% figure(7)
% subplot(2,1,1);
% plot(tn,abs(MZI_bar).^2);
% xlim([10,15]);
% xlabel('time[ns]');
% ylabel('Intensity[mW]');
% % saveas(gcf,'15.jpg');
% 
% % ylim([-4,4]);
% subplot(2,1,2);
% plot(tn,unwrap(angle(MZI_bar)));
% ylim([-4,4]);
% xlim([10,15]);
% xlabel('time[ns]');
% ylabel('Phase[rad]');

figure(30)
subplot(2,1,1)
plot(tn-5.4,abs(MZI_bar).^2 * 10^3);
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Intensity[mW]');
subplot(2,1,2)
plot(tn-5.4,angledef(MZI_bar));
xlim([0,0.5]);
xlabel('time[ns]');
ylabel('Phase[rad]');

%Q factor

MZI_1 =[];
MZI_0 =[];

%Q??�l??�̂�??�߂̑S??�|??�C??�??�??�g??�??�??�o
% for ii = 1:1:size(t,2)
%     if MZI_cross(1,ii) < 0
%         MZI_1 = [MZI_1;MZI_cross(1,ii)];
%     end
%     if MZI_cross(1,ii) > 0
%          MZI_0 = [MZI_0;MZI_cross(1,ii)];
%     end
%     
% end

point_div = 10;

%%%%%%%%??�ʑ�??�L??�??�??�??�̕�??�U??�̌v??�Z%%%%%%%%%%%%%%%%%%%%%%
phase_phai_1 = [];
phase_phai_0 = [];
capture_sig = [];

for i_c = 1:1:size(MZI_bar,1)/bit_d
    capture_sig(i_c,1) = MZI_bar(i_c * bit_d - bit_d/2);
end


for i = 1:1:size(capture_sig,1)
        if abs(angle(capture_sig(i))) > 1/2 * pi()
%         if(angle(capture_sig(i))) > 0
             phase_phai_1 =[phase_phai_1;capture_sig(i)];
        else
             phase_phai_0 =[phase_phai_0;capture_sig(i)];
%         end
%     for ii = -point_div:1:point_div
        end
%     end
end

sigma_phai_out_zero = sum((0-abs(angle(phase_phai_0))).^2)/size(phase_phai_0,1);
sigma_phai_out_one = sum((pi()-abs(angle(phase_phai_1))).^2)/size(phase_phai_1,1);

sigma_phai_out = (sigma_phai_out_one + sigma_phai_out_zero)/2;
%%%%%%%%??�ʑ�??�L??�??�??�??�̕�??�U??�̌v??�Z%%%%%%%%%%%%%%%%%%%%%%



% sigma_phai_out = var(angle(phase_phai) + pi());

% %??�s??�[??�N??�_??�܂�??�̃|??�C??�??�??�g??�??�??�o
% 
% point_div = 70;
% 
% for i = 1:1:patt_len
%     for ii = -point_div:1:point_div
%         if MZI_bar(i * bit_d/2 + ii,1) < 0
%         MZI_1 = [MZI_1;MZI_bar(1,i * bit_d/2 + ii)];
%         end
%     end
%     if MZI_bar(1,i * bit_d/2 + ii) > 0
%         MZI_0 = [MZI_0;MZI_bar(1,i * bit_d/2 + ii)];
%     end
% end
%     
% MZI_mean_0 = mean(abs(MZI_0).^2);
% MZI_var_0  = var(abs(MZI_0).^2);
% MZI_mean_1 = mean(abs(MZI_1).^2);
% MZI_var_1  = var(abs(MZI_1).^2);
% 
% patt_len;
% bit_d;
% 
% Q = 10* log10((abs(MZI_mean_1 - MZI_mean_0))/(MZI_var_1+MZI_var_0))
% disp('[dB]')

%??�R??�??�??�X??�^??�??�??�[??�V??�??�??�??�f
% MZI_bar = signal_gain(MZI_bar,0);
conste_temp = MZI_bar;
% Generate_Constellation(conste_temp,bit_d);
% Generate_Eyepattern(conste_temp,t,bit_d,patt_len,dat_len);
clear conste_temp


 
[E_out,Q] = DPSK_receiver(MZI_bar,bit_d);
% Generate_Eyepattern(E_out,bit_d);
Q
sigma_phai_out
% load('make_colormap.mat','Q_map','sigma_map');
Q_map =[Q_map;[i_alpha,i_i_s*1e3,Q]];
sigma_map = [sigma_map;[i_alpha,i_i_s*1e3,sigma_phai_out]];
% save('make_colormap.mat','Q_map','sigma_map');

% end         %sa
% end         %sa

% hMod = modem.pskmod('M', 4, 'PhaseOffset', pi/4);
% hScope = commscope.ScatterPlot
% hScope.Constellation = hMod.Constellation;
% hScope.SamplesPerSymbol = 160;
% update(hScope, MZI_bar)
% scatter(abs(MZI_1).^2 .* cos(angle(MZI_1)),abs(MZI_1).^2 .* sin(angle(MZI_1)));
% scatter(abs(MZI_0).^2 .* cos(angle(MZI_0)),abs(MZI_0).^2 .* sin(angle(MZI_0)));
% scatter(cos(angle(MZI_bar)),sin(angle(MZI_bar)));
% 
% 
% xlim([-1,1]);
% ylim([-1,1]);
% 
% [x,y,z] = cylinder(1,100);
% z(:) = 0; % ??�??�??�??�??�f??�[??�^??�??�0??�ɐݒ�
% surf(x,y,z)
% view(2)
% axis([-1,1,-1,1])
% axis square
% grid off  
% xlabel('I-Phase');
% ylabel('Q-Phase');
% hold off

% xlabel('t (second)')
% ylabel('Intensity')
% saveas(gcf,'E_in.jpg')

%
%%H_fiber = exp(-1i *(0.5*Dw*f.^2));
%
%%LCOS??�̉𑜓x??�ɂ�鐧??�??�
%Lcos_delta=4;
%for index = 1:Lcos_delta:dat_size_P
%        for roop = 1:1:Lcos_delta-2
%                H_fiber(index+roop,1)=H_fiber(index,1);
%        end
%end
%??�T??�C??�Y??�??�??�??�
%if size(H_fiber,1)-size(f,1) ~= 0
%        H_fiber=H_fiber(1:size(f),1);
%end
%
%LCOS??�̊K??�??�??�??�??�ɂ�鐧??�??�
%Lcos_power_tmp=abs(H_fiber);
%Lcos_phase_tmp=angle(H_fiber);
%for i =1:1:dat_size_P
%        for ii = 1:1:Lcos_grad+1
            
%                if Lcos_phase_tmp(i,1)>0
%                    if Lcos_phase_tmp(i,1) < 2*pi/Lcos_grad*ii
%                        Lcos_phase_tmp(i,1) = 2*pi/Lcos_grad*(ii-1);
%                       break;
%                    end
%                else
%                    if Lcos_phase_tmp(i,1) > -2*pi/Lcos_grad*ii
%                        Lcos_phase_tmp(i,1) = -2*pi/Lcos_grad*(ii-1);
%                        break;
%        
%                    end
%                end
%        end
%end          



%H_fiber =exp(1i*Lcos_phase_tmp);

%??�T??�C??�Y??�??�??�??�
%if size(H_fiber,1)-size(f,1) ~= 0
%        H_fiber=H_fiber(1:size(f),1);
%end
%size(H_fiber)
%size(f)

%plot(f,angle(H_fiber));
%xlim([0,3*10^10]);
%ylim([-pi,pi]);

%F_out = F_in .* H_fiber;

%E_out = ifft(F_out);

%complex_waveform = E_out;

% %??�f??�[??�^??�擾??�p
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(t,E_out)
% xlabel('t (second)')
% ylabel('Intensity')
% saveas(gcf,'E_out.jpg')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_out = MZI_bar;
Power = (abs(E_out)).^2;

Phase    = (angle(E_out));
% disp('--------------------------------------------------------------------')
%Phase    = unwrap(angle(E_out) + pi*0.5*((-1).^(n+1)+1)');

% Phase    = (angle(E_out));

