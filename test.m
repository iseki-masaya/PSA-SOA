clear;                               %initialize
close all;

disp('     ===================================     ')
disp('     =====    Dual-pump SOA PSA    =====     ')
disp('     =====      Trans. & Regene.   =====     ')
disp('     =====        ver. 1.10        =====     ')
disp('     ===================================     ')

tic

parameter;                           %setting of parameter (constants)

bit_d           = 10;                       %1?r?b?g?f?[?^????????
del_t           = 1 / (speed * bit_d);      %????????[s]

input_condition;                     %setting of input conditions

DPSK_waveform_generation;            %making optical waveform

Fiber_characteristics;              %setting of fiber characteristics

E_out = sqrt(P_input_sig).*exp(1i * Phase_input_sig);

E_out = awgn(E_input_sig,15,'measured');

Power_out = abs(E_input_sig).^2;
Phase_out = angle(E_input_sig);

E_out = E_input_sig;
F_out = fft(E_out);
df = 1/del_t;
f = 0:df/size(E_out,1):df;
f = f';
f(end,:) = [];
%plot(f,abs(F_out).^2/max(abs(F_out).^2))
xlim([0,df/2]);

%[inp_sig,Q_in]=DPSK_receiver(E_out,bit_d);

tic
n = 10;

for i = 1:1:n
	%%%%%%%%%%%%%%%%%%%%%%%%%SPM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[E_out_SPM, Power_out_SPM, Phase_out_SPM] =  SPM(Power_out,del_t,...
	    Phase_out,D_SSMF,S_SSMF,100e3,1565e-9, 1550e-9); 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	ASE = ASE_noise(2);
	%E_ASE = randn(size(E_out_SPM,1),1).*sqrt(ASE);
	E_ASE = sqrt(ASE);
    E_ASE_phase = rand(size(E_out_SPM,1),1)*2*pi();

	E_out_SPM = E_out_SPM + E_ASE.*exp(Im*E_ASE_phase);
	
	% SNR = 10*log10(mean(abs(E_out_SPM)).^2/ASE);
	% E_out_SPM = awgn(E_out_SPM,SNR,'measured');
	Power_out_SPM = abs(E_out_SPM).^2;
	Phase_out_SPM = angle(E_out_SPM);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%???U??%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[E_out, Power_out,  Phase_out] =  Fiber_Transmission(Power_out_SPM,...
 	   del_t, Phase_out_SPM,D_DCF,S_DCF,100e3,1565e-9, 1550e-9); 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	toc
	t = get_time(size(E_out))*10^9;

end

%Ideal filter

%E_out_filtered = fft(E_out);
%E_out_filtered(275:5200,1) = 0;
%E_out_filtered = ifft(E_out_filtered);

E_out_filtered = E_out;

disp('degraded signal')

Generate_Constellation(E_out_filtered,bit_d);
line = findobj('type','line');
line_obj_degenerated=get(line(3));
save line_obj_degraded.mat line_obj_degenerated;
%plot(line_obj_degenerated.XData,line_obj_degenerated.YData);

%EVM measurement
size_data=numel(line_obj_degenerated.XData);
error_vector_deg_0=0;
error_vector_deg_1=0;
mean_0_x = 0;
mean_0_y = 0;
mean_1_x = 0;
mean_1_y = 0;
count0=0;
count1=0;

for iiii=1:1:size_data
    if line_obj_degenerated.XData(1,iiii) > 0
        count0 = count0 + 1;
    else
        count1 = count1 + 1;
    end
end

vector0=zeros(1,count0);
vector1=zeros(1,count1);

count0=0;
count1=0;

for iiii=1:1:size_data
    if line_obj_degenerated.XData(1,iiii) > 0
        mean_0_x = mean_0_x + line_obj_degenerated.XData(1,iiii);
        mean_0_y = mean_0_y + line_obj_degenerated.YData(1,iiii);
        count0 = count0 + 1;
        vector0(1,count0) = line_obj_degenerated.XData(1,iiii) + Im * line_obj_degenerated.YData(1,iiii);
    else
        mean_1_x = mean_1_x + line_obj_degenerated.XData(1,iiii);
        mean_1_y = mean_1_y + line_obj_degenerated.YData(1,iiii);
        count1 = count1 + 1;
        vector1(1,count1) = line_obj_degenerated.XData(1,iiii) + Im * line_obj_degenerated.YData(1,iiii);
    end
end

ref_vector_deg_0 = (mean_0_x + Im * mean_0_y) / count0;
ref_vector_deg_1 = (mean_1_x + Im * mean_1_y) / count1;

count0=0;
count1=0;

for iiii=1:1:size_data
    if line_obj_degenerated.XData(1,iiii) > 0
        error_vector_deg_0 = error_vector_deg_0 + (abs(line_obj_degenerated.XData(1,iiii)+...
            Im*line_obj_degenerated.YData(1,iiii)-ref_vector_deg_0))^2;
        count0=count0+1;
    else
        error_vector_deg_1 = error_vector_deg_1 + (abs(line_obj_degenerated.XData(1,iiii)+...
            Im*line_obj_degenerated.YData(1,iiii)-ref_vector_deg_1))^2;
        count1=count1+1;
    end
end
error_vector_deg_0 = sqrt(error_vector_deg_0/count0)/abs(ref_vector_deg_0);
error_vector_deg_1 = sqrt(error_vector_deg_1/count1)/abs(ref_vector_deg_1);
error_vector_deg=(error_vector_deg_0 + error_vector_deg_1)/2

phase_offset_degraded_sig = abs((angle(ref_vector_deg_1)+(angle(ref_vector_deg_0)-pi))/2);

%ref_vector_deg_0 = mean_0 * exp(Im * phase_offset_degraded_sig);
%ref_vector_deg_1 = mean_1 * exp(Im * phase_offset_degraded_sig);

E_out_filtered=E_out_filtered * exp(Im * phase_offset_degraded_sig);

%Bit error rate mesurement
Q_phase = 10*log10(0.87 * (pi() / (std(unwrap(angle(vector0))) + std(unwrap(angle(vector1))))))
%Q_phase = 20 * log10(Q_phase)
BER = erfc(Q_phase/sqrt(2))

Generate_Constellation(E_out_filtered,bit_d);
%Generate_Eyepattern(abs(E_out).^2,bit_d);
Generate_Eyepattern((E_out_filtered),bit_d);

[output_sig,Q]=DPSK_receiver(E_out_filtered,bit_d);
sigma = get_sigma(E_out_filtered,bit_d);
%Q

toc
