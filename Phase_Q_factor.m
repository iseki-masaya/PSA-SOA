function [ phase_q ] = Phase_Q_factor(E_in,bit_d )

phase_phai_1 = [];
phase_phai_0 = [];
capture_sig = [];

for i_c = 1:1:size(E_in,1)/bit_d
    capture_sig(i_c,1) = angle(E_in(i_c * bit_d - bit_d/2));
end

capture_sig=unwrap(capture_sig);

for i = 1:1:size(capture_sig,1)
        if abs((capture_sig(i))) > 1/2 * pi()
%         if(angle(capture_sig(i))) > 0
             phase_phai_1 =[phase_phai_1;capture_sig(i)];
        else
             phase_phai_0 =[phase_phai_0;capture_sig(i)];
%         end
%     for ii = -point_div:1:point_div
        end
%     end
end

%sigma_phai_out_zero = sum((0-abs(angle(phase_phai_0))).^2)/size(phase_phai_0,1);
%sigma_phai_out_one = sum((pi()-abs(angle(phase_phai_1))).^2)/size(phase_phai_1,1);

%sigma_phai_out = (sigma_phai_out_one + sigma_phai_out_zero)/2;
%GET_SIGMA Summary of this function goes here
%   Detailed explanation goes here


mean_0 = mean(phase_phai_0);
var_0  = sqrt(var(phase_phai_0));
mean_1 = mean(phase_phai_1);
var_1  = sqrt(var(phase_phai_1));

phase_Q = 20* log10((abs(mean_1 - mean_0))/(var_1+var_0));
phase_Q_factor = [phase_Q_factor;phase_Q];
phase_q = max(phase_Q_factor);

end