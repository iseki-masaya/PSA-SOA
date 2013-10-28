clear;                               %initialize

parameter;                           %setting of parameter (constants)

input_condition;                     %setting of input conditions

DPSK_waveform_generation;            %making optical waveform

%Fiber_characteristics;              %setting of fiber characteristics

E_out = sqrt(P_input_sig).*exp(1i * Phase_input_sig);
E_out = awgn(E_input_sig,7,'measured');

Generate_Constellation(E_out,bit_d);
Generate_Eyepattern(abs(E_out).^2,bit_d);

[output_sig,Q]=DPSK_receiver(E_out,bit_d);
sigma = get_sigma(E_out,bit_d)
Q

%P_input_sig = abs(E_out).^2;
%Phase_input_sig = Phase_out;
