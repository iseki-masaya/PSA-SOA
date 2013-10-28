%function [output_sig ] = SESAM(input_sig)
function [output_sig] = SESAM(input_sig,i_f,i_b,i_alpha,i_i_s)

input_phase = unwrap(angle(input_sig));
% input_power = abs(input_sig).^2 / (r^2 * pi());
input_power = abs(input_sig).^2;
% r_1 = 0.7;
% r_2 = 0.945;
r_1 = i_f;
r_2 = i_b;

R_1 = r_1;
R_2 = r_2;
n_eff_k_0_l  = 2 * pi();
% n_eff_k_0_l  = pi()/2;

alpha_0 = i_alpha;
 % alpha_0 = 0.053/2;   
%         I_s = 0.1* 10^(-4);
I_s = i_i_s;
          
% I_change = ((14 * alpha_0 * I_s)^4-9 * I_s)/9;
alpha_l =alpha_0 ./(1+input_power/I_s);
T_a = exp(-2*alpha_l);

R_b_eff = T_a * R_2;

I_t = (1-R_1)*(1-R_2)*T_a./((1-sqrt(R_1*R_2.*T_a)).^2 + 4 * sqrt(R_1*R_2.*T_a)* sin(n_eff_k_0_l)^2);

I_r = ((sqrt(R_1)-sqrt(R_2 .* T_a)).^2 + 4*sqrt(R_1*R_2.*T_a)*sin(n_eff_k_0_l))./((1-sqrt(R_1*R_2.*T_a)).^2+4*sqrt(R_1*R_2.*T_a)*sin(n_eff_k_0_l).^2);

% output_sig = sqrt(input_power .* I_r* (r^2 * pi())) .* exp(1i.*input_phase); 
output_sig = sqrt(input_power .* I_r) .* exp(1i.*input_phase); 
figure(15);
subplot(2,1,1)
    plot(R_b_eff,I_r,'r');
    xlabel('R^b_{eff} (function of input power)');
    ylabel('Cavity reflectivity at resonance');
subplot(2,1,2)
    plot(input_power,R_b_eff,'b');
    xlabel('input power[mW]');
    ylabel('R^b_{eff} (function of input power)');
    
%     
% % T_dif = (max(T_a)-exp(-2*7*alpha_0/2)) * 100
% % 
% % if T_max < T_dif
% %     T_max = T_dif;
% %     eval(['print -djpeg -r70 ./temp_cal/alpha-' num2str(i) '.jpeg']);
% % 
% % end
% 
% % end
% 
% 
% 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%?O???t???p?i???o??j%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % r_1 = 0.9;
% % % r_2 = 0.945;
% r_1 = i_f;
% r_2 = i_b;
% R_1 = r_1;
% R_2 = r_2;
% 
% I_r = ((sqrt(R_1)-sqrt(R_2 .* T_a)).^2 + 4*sqrt(R_1*R_2.*T_a)*sin(n_eff_k_0_l))./((1-sqrt(R_1*R_2.*T_a)).^2+4*sqrt(R_1*R_2.*T_a)*sin(n_eff_k_0_l).^2);
% I_t = (1-R_1)*(1-R_2)*T_a./((1-sqrt(R_1*R_2.*T_a)).^2 + 4 * sqrt(R_1*R_2.*T_a)* sin(n_eff_k_0_l)^2);
% 
% % output_sig = sqrt(input_power).*exp(1i.*input_phase);
% hold on
% plot(R_b_eff,I_t);
% 
% eval(['print -djpeg -r70 ./rev_jpeg/fornt-' num2str(i_f) 'back-' num2str(i_b) '.jpeg']);
% hold off
end

