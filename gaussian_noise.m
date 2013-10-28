function [n1,n2,n3,n4,n5,n6] = gaussian_noise()

parameter;

c               = 2.99792458e8;
term = 300;                       %����
eta_i           = 0.8;      %�O���ʎq��
n_sp            = 1.3;        %���R��o�W��
m_t             = 2;        %���[�h��
lambda          = 1550*10^(-9); % �g��
B_e             = 10*10^9;                        % �M���ш�
%dat_len;
Res             = 50;        %�O����R[��]
B_opt           = 100*10^9;  %BPF�̑ш�
q_unit          = 1.60217646e-19; 
h               = 6.6260693e-34; 
boltzmann       = 1.3806523e-23; 
R = 50;
input_level_sig = -10;
P_s             = 10^(input_level_sig / 10) * 10^(-3);
bit_d           = 160;
n = 1;
G_0 = 10;
G_1 = 100;   %(20dB)
SIG =G_0; 
ASE =SIG-1+(G_1-1)*n;

RIN = -135-10 * log10((P_s/10^(-3))^1);     %[dB]

N_shot_p = 2 * q_unit * (q_unit * eta_i/(h * c/lambda)) * P_s * B_e * SIG;    %���x
N_spshot =  2 * q_unit * (q_unit * eta_i * B_e * n_sp * m_t * B_opt) * ASE;
N_th_p = 4 * boltzmann * term / Res * B_e/R;                              %���x
N_ssp = 4 * (q_unit * eta_i * P_s/(h * c/lambda))* (q_unit * eta_i * n_sp * B_e)*SIG*ASE;
N_spsp = 4 * (q_unit * eta_i * n_sp)^2 * m_t * B_opt * B_e * ASE^2;
N_RIN = sqrt(10^(RIN /10*B_e)*(P_s/2)^2) * SIG;

% n2 = 4
% n1 = n2-1;
% G1 = 100*n1;
% G2 = 100*n2;
G = 100;

N_ASE1 = h * c/lambda * (G-1) * B_e;

i_s =  (q_unit * eta_i/(h * c/lambda*G))^2;
i_s =  (q_unit * eta_i/(h * c/lambda*G))^2;

i_b =  2 * (q_unit * eta_i/(h * c/lambda))^2*G*P_s*N_ASE1;

i_sh1 = 2 * q_unit * (q_unit * eta_i/(h * c/lambda)) * P_s * B_e;

SNR_out = 10*log10(i_s/(i_sh1+i_b));
SNR_in = 10*log10(i_s/(i_sh1));
NF = SNR_in/SNR_out;

n1 = N_ssp;
n2 = N_spsp;
n3 = N_RIN;
n4 = N_th_p ;
n5 = N_spshot;
n6 = N_shot_p;
% I = (q_unit * eta_i/(h * c/lambda)) * P_s;
% 
% 
% 
% P_s
% I
% R_L_I = Res * I
% N_2eI = 2 * q_unit * I
% N_4kt = 4 * boltzmann * term / Res
% SNR = (I)^2/(N_shot_p/2 + N_th_p/2)
% SNR_dB = 10 * log10(SNR)
% 
% %�z���C�g�m�C�Y�����肵���M�G��
% ran = rand(num,1) * sqrt(N_th_p);
% ranf = fft(ran);
% ranf_i = abs(ranf) * N_th_p;
% ranf_p = angle(ranf);
% N_th = ifft(ranf_i .* exp(1i * ranf_p));
% 
% %�|�A�\�����z�����肵���V���b�g�G���i�v�C��)
% 
%  ran = rand(num,1) * sqrt(N_shot_p);
% ranf = fft(ran);
% ranf_i = abs(ranf) * N_shot_p;
% ranf_p = angle(ranf);
% N_shot = ifft(ranf_i .* exp(1i * ranf_p));
% 
% %J = N_th + N_shot;
% J = SNR_dB;

