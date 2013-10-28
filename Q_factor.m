
parameter;                           %setting of parameter (constants)

input_condition;                     %setting of input conditions

DPSK_waveform_generation;            %making optical waveform

Fiber_characteristics;               %setting of fiber characteristics

% E_out = signal_gain(E_out,-10);

% F_out = fft(E_out);
% df = 1/del_t;
% f = 0:df/size(E_out,1):df;
% f = f';
% f(end,:) = [];
% plot(f,abs(F_out).^2/max(abs(F_out).^2))
% xlim([0,df/2]);
Generate_Constellation(E_out,bit_d);
Generate_Eyepattern(abs(E_out).^2,bit_d);
% 
[output_sig,Q]=DPSK_receiver(E_out,bit_d);
sigma = get_sigma(E_out,bit_d)
Q
% output_sig = output_sig*150;
% Generate_Eyepattern(output_sig,160);
% 
% output_sig = retiming(output_sig,160);
% figure(1)
% x = -0.5:0.01:0.5;
% hist(output_sig,x);
% 
% 
%{
 P_input_sig = abs(E_out).^2;
 Phase_input_sig = Phase_out;
 
figure(3)
subplot(2,1,1)
plot(t-1,abs(E_out).^2*1e3);
xlabel('time[ns]');
ylabel('Intensity[mW]');
xlim([0,0.5])

subplot(2,1,2)
plot(t-1,Phase_out);
xlabel('Phase[rad]');
ylabel('Intensity[mW]');
xlim([0,0.5])
%}
% Generate_Constellation(E_out,bit_d);
% Generate_Eyepattern(abs(E_out).^2,bit_d);
% 
% [output_sig,Q]=DPSK_receiver(E_out,bit_d);
% 
% output_sig = output_sig*150;
% Generate_Eyepattern(output_sig,bit_d);
% 
% output_sig = retiming(output_sig,bit_d);
% figure(1)
% x = -0.5:0.01:0.5;
% hist(output_sig,x);
% 
% 
