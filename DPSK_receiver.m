function [output_sig,Q_factor] = DPSK_receiver(input_sig,bit_d)

parameter
R = 500;		%?½t?½?½?½?½?½g?½G?½?½?½h?½?½R
eta_i = 0.8;        %?½O?½?½?½ÊŽq?½?½
lambda          = 1550*10^(-9); % ?½g?½?½
B_e             = 10*10^9;    
q_unit          = 1.60217646e-19;           %?d???f??

% figure(60)
% plot(real(input_sig));

input_sig1 = [input_sig;zeros(bit_d,1)]/2;   %top arm separating original singnal 
input_sig2 = [zeros(bit_d,1);input_sig,]/2;     %bottom arm separating original singaland operate a 1bit delay

temp_sig1 = (input_sig1 - input_sig2);  %interference between top and bottom arms
temp_sig2 = (input_sig1 + input_sig2);

 for i = 1:1:bit_d
    temp_sig1(1,:) = [];
    temp_sig1(end,:) = [];
    temp_sig2(1,:) = [];
    temp_sig2(end,:) = [];
 end
 
recieve_sig = q_unit * eta_i/(h * c/lambda )* (abs(temp_sig1).^2 - abs(temp_sig2).^2);

dat_len = size(temp_sig1,1);
patt_len = dat_len/bit_d;
noise_size= size(temp_sig1,1);

[~,~,~,N_th,~,~] = gaussian_noise();

I_th = (randn(noise_size,1) .* sqrt(N_th));      %?½G?½?½?½d?½?½

N_cir = sqrt(B_e) /100 * B_e * q_unit;
I_cir = (randn(noise_size,1) .* N_cir);
I_shot = shot_noise(recieve_sig);
noise_t = I_th + I_shot;

output_sig = recieve_sig + noise_t;

%Generate_Eyepattern((output_sig),bit_d);

% %%%%%%%?½?½?½?½M?½t?½?½?½?½?½g?½G?½?½?½h%%%%%%%%%%%%
% gm = 4*10e-1;
% re = 1*10e2;
% RL = 50;
% gds = 1*10e-13;
% R = 500;
% Cc = 2*10e-13;
% Ct = 1*10e-13;
% Cds = 1*10e-13;
% gamma = 0.03;		%?½G?½?½?½áŒ¸?½W?½?½
% boltzmann       = 1.3806523e-23;
% Tn = 274;
% Bn = 10*10e9;
% fc = 2.99792458e8/lambda;
% omega = 2 * pi() * fc * 10^9;
% 
% Y = 1/R + 1i * omega * Ct + (1i * omega * Cc)/(re*(1i * omega * Cc) + 1);
% 
% 
% FI = fft(recieve_sig);
% 
% 
% E0 = FI /Y * (1i * omega * Cc)/(re * (1i * omega * Cc)+1);
% 
% 
% pfet = noise_th*sqrt(10^3);
% pFet = fft(noise_th);
% 
% 
% E0n = pFet./((re + 1/(1/R + 1i * omega * Ct))*(1i * omega * Cc)+1);
% 
% F_out = gm * (E0 + E0n)./(1/RL + gds + 1i * omega * Cds);
% output_sig = ifft(F_out);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num = size(output_sig,1);

sig_1_level = [];
sig_0_level = [];
q_cal =output_sig ;

q_cal = q_cal/mean(abs(q_cal));

Q_factor =[];

for iii = 1:1:(bit_d/2-1)
    point_div = iii;
    % point_div = 10;


    for i = 1:1:num/bit_d
        for ii = -point_div:1:point_div
            if q_cal(i * bit_d -bit_d/2 + ii,1) > 0
                sig_1_level = [sig_1_level;q_cal(i * bit_d -bit_d/2 + ii,1)];
            else
                sig_0_level = [sig_0_level;q_cal(i * bit_d -bit_d/2 + ii,1)];
            end
        end
    end


    % for i = 1:1:num
    %     if q_cal(i,1) > 0
    %         sig_1_level = [sig_1_level;q_cal(i,1)];
    %     else
    %         sig_0_level = [sig_0_level;q_cal(i,1)];
    %     end
    % end

    mean_0 = mean(sig_0_level);
    var_0  = sqrt(var(sig_0_level));
    mean_1 = mean(sig_1_level);
    var_1  = sqrt(var(sig_1_level));

    Q = 20* log10((abs(mean_1 - mean_0))/(var_1+var_0));
    Q_factor = [Q_factor;Q];
    Q_factor = max(Q_factor);
    % figure(70)
    % plot(Q_factor)
end

