function [ E_out ] = signal_amp( E_in,dB)

x = 10^(dB/10);
P_in = abs(E_in).^2;
P_in_phase = angle(E_in);

P_out = P_in*x;
E_out = sqrt(P_out) .* exp(1i * P_in_phase);
%SIGNAL_AMP Summary of this function goes here
%   Detailed explanation goes here


end

