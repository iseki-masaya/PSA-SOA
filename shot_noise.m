function [I_shot] = shot_noise(I_in)
%SHOT_NOISE この関数の概要をここに記述
parameter;

c               = 2.99792458e8;
term = 300;                       %室温
eta_i           = 0.8;      %外部量子効率
n_sp            = 1.3;        %自然放出係数
m_t             = 2;        %モード数
lambda          = 1550*10^(-9); % 波長
B_e             = 1*10^9;                        % 信号帯域
Res             = 50;        %外部抵抗[Ω]
B_opt           = 100*10^9;  %BPFの帯域
q_unit          = 1.60217646e-19; 
h               = 6.6260693e-34; 
boltzmann       = 1.3806523e-23;   
noise_size = size(I_in,1);
I_mean = mean(abs(I_in));

N_shot_p = 2 * q_unit * I_mean * B_e;    %今度
I_shot = randn(noise_size,1) .* sqrt(N_shot_p);   
end

