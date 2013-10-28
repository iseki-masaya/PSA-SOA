function [I_shot] = shot_noise(I_in)
%SHOT_NOISE ���̊֐��̊T�v�������ɋL�q
parameter;

c               = 2.99792458e8;
term = 300;                       %����
eta_i           = 0.8;      %�O���ʎq����
n_sp            = 1.3;        %���R���o�W��
m_t             = 2;        %���[�h��
lambda          = 1550*10^(-9); % �g��
B_e             = 1*10^9;                        % �M���ш�
Res             = 50;        %�O����R[��]
B_opt           = 100*10^9;  %BPF�̑ш�
q_unit          = 1.60217646e-19; 
h               = 6.6260693e-34; 
boltzmann       = 1.3806523e-23;   
noise_size = size(I_in,1);
I_mean = mean(abs(I_in));

N_shot_p = 2 * q_unit * I_mean * B_e;    %���x
I_shot = randn(noise_size,1) .* sqrt(N_shot_p);   
end

