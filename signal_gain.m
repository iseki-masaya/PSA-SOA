function [ E_out ] = signal_gain( E_in,dBm)

x_mW = 10^(dBm/10);
P_in = abs(E_in).^2;
P_in_phase = angle(E_in);

P_out = P_in/max(P_in) * x_mW * 1e-3;
E_out = sqrt(P_out) .* exp(1i * P_in_phase);

% E_out = E_in/max(E_in) * sqrt(x_mW * 1e-3);

%SIGNAL_GAIN ���̊֐��̊T�v�������ɋL�q
%   �ڍא����������ɋL�q


end

