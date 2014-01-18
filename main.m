clear;                               %initialize
close all;

disp('     ===================================     ')
disp('     =====    Dual-pump SOA PSA    =====     ')
disp('     =====      Trans. & Regene.   =====     ')
disp('     =====        ver. 1.10        =====     ')
disp('     ===================================     ')

tic

parameter;                           %setting of parameter (constants)

input_condition;                     %setting of input conditions

DPSK_waveform_generation;            %making optical waveform

Fiber_characteristics;              %setting of fiber characteristics

Power_out = abs(E_input_sig).^2;
Phase_out = angle(E_input_sig);

E_out = E_input_sig;
F_out = fft(E_out);
df = 1/del_t;
f = 0:df/size(E_out,1):df;
f = f';
f(end,:) = [];
plot(f,abs(F_out).^2/max(abs(F_out).^2))
xlim([0,df/2]);
toc

tic
n = 10;

for i = 1:1:n
	%%%%%%%%%%%%%%%%%%%%%%%%%SPM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[E_out_SPM, Power_out_SPM, Phase_out_SPM] =  SPM(Power_out,del_t,...
	    Phase_out,D_SSMF,S_SSMF,100e3,1565e-9, 1550e-9); 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	ASE = ASE_noise(2);
	E_ASE = sqrt(ASE);
    E_ASE_phase = rand(size(E_out_SPM,1),1)*2*pi();

	E_out_SPM = E_out_SPM + E_ASE.*exp(Im*E_ASE_phase);
	
	Power_out_SPM = abs(E_out_SPM).^2;
	Phase_out_SPM = angle(E_out_SPM);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%???U??%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[E_out, Power_out,  Phase_out] =  Fiber_Transmission(Power_out_SPM,...
 	   del_t, Phase_out_SPM,D_DCF,S_DCF,100e3,1565e-9, 1550e-9); 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	toc
    disp(['count : ',i]);
	t = get_time(size(E_out))*10^9;
end

E_out_filtered = E_out;

disp('degraded signal')

Generate_Constellation(E_out_filtered,bit_d);
line = findobj('type','line');
line_obj_degenerated=get(line(3));
save line_obj_degraded.mat line_obj_degenerated;
plot(line_obj_degenerated.XData,line_obj_degenerated.YData);

%EVM measurement
size_data = numel(line_obj_degenerated.XData);
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
error_vector_deg = (error_vector_deg_0 + error_vector_deg_1)/2

phase_offset_degraded_sig = abs((angle(ref_vector_deg_1)+(angle(ref_vector_deg_0)-pi))/2);

E_out_filtered=E_out_filtered * exp(Im * phase_offset_degraded_sig);

%Generate_Constellation(output_sig,bit_d);
%Generate_Eyepattern(output_sig,bit_d);

%P_input_sig = abs(E_out).^2;
%Phase_input_sig = Phase_out;

%Preamplification morep precise concideration is required
preamp = 10;                       %gain [dB]
E_out_filtered = 10^(preamp/10)*E_out_filtered;

B = 100 * 10^9;                     %ASE band width
ase_at_preamp = (h * c/lam_sig)*(10^(preamp/10)-1)*B;  %Assuming nsp as 1

input_sig = (abs(E_out_filtered)).^2.*exp(Im.*angle(E_out_filtered))+ase_at_preamp.*exp(Im*randn());

TMM                                 %Transfer Matrix Method

disp('regenerated signal')

num = size(A_sig,1);

E_out_sig = zeros(num-div_n,1);

for i = 1:1:num-div_n
    E_out_sig(i,1)=A_sig(i+div_n,div_n);
end

E_out_sig_filtered = E_out_sig;

Generate_Constellation(E_out_sig_filtered,bit_d);

line = findobj('type','line');
line_obj_regenerated=get(line(5));
plot(line_obj_degenerated.XData,line_obj_degenerated.YData);

size_data=numel(line_obj_regenerated.XData);
error_vector_reg_0=0;
error_vector_reg_1=0;
mean_0_x = 0;
mean_0_y = 0;
mean_1_x = 0;
mean_1_y = 0;
count0=0;
count1=0;

for iiii=1:1:size_data
    if line_obj_regenerated.XData(1,iiii) > 0
        count0 = count0 + 1;
    else
        count1 = count1 + 1;
    end
end

vector0_reg=zeros(1,count0);
vector1_reg=zeros(1,count1);

count0=0;
count1=0;

for iiii=1:1:size_data
    if line_obj_regenerated.XData(1,iiii) > 0
        mean_0_x = mean_0_x + line_obj_regenerated.XData(1,iiii);
        mean_0_y = mean_0_y + line_obj_regenerated.YData(1,iiii);
        count0 = count0 + 1;
        vector0_reg(1,count0) = line_obj_regenerated.XData(1,iiii) + Im * line_obj_regenerated.YData(1,iiii);
    else
        mean_1_x = mean_1_x + line_obj_regenerated.XData(1,iiii);
        mean_1_y = mean_1_y + line_obj_regenerated.YData(1,iiii);
        count1 = count1 + 1;
        vector1_reg(1,count1) = line_obj_regenerated.XData(1,iiii) + Im * line_obj_regenerated.YData(1,iiii);
    end
end

ref_vector_reg_0 = (mean_0_x + Im * mean_0_y) / count0;
ref_vector_reg_1 = (mean_1_x + Im * mean_1_y) / count1;

phase_offset_regenerated_sig = abs((angle(ref_vector_reg_1)+(angle(ref_vector_reg_0)-pi))/2);

count0=0;
count1=0;

for iiii=1:1:size_data
    if line_obj_regenerated.XData(1,iiii) > 0
        error_vector_reg_0 = error_vector_reg_0 + (abs(line_obj_regenerated.XData(1,iiii)+...
            Im*line_obj_regenerated.YData(1,iiii)-ref_vector_reg_0))^2;
        count0=count0+1;
    else
        error_vector_reg_1 = error_vector_reg_1 + (abs(line_obj_regenerated.XData(1,iiii)+...
            Im*line_obj_regenerated.YData(1,iiii)-ref_vector_reg_1))^2;
        count1=count1+1;
    end
end
error_vector_reg_0 = sqrt(error_vector_reg_0/count0)/abs(ref_vector_reg_0);
error_vector_reg_1 = sqrt(error_vector_reg_1/count1)/abs(ref_vector_reg_1);
error_vector_reg=(error_vector_reg_0 + error_vector_reg_1)/2

Generate_Eyepattern((E_out_sig_filtered),bit_d);
[output_sig,Q]=DPSK_receiver(E_out_sig_filtered,bit_d);
sigma = get_sigma(E_out_sig_filtered,bit_d);

ASE_spectral = zeros(30,1);
ASE_spectral_axis = zeros(30,1);
for ii=1:30
    ASE_spectral(ii,1)=10 * log10(ASE_spectral(ii,1)/1e-3);
    ASE_spectral_axis(ii,1)=(1503+(ii-1)*4)*1e-9;
end

%signals spectral
fft_sig = fftshift(fft(P_for_sig(:,div_n)));
fft_pu1 = fftshift(fft(P_for_pu1(:,div_n)));
fft_pu2 = fftshift(fft(P_for_pu2(:,div_n)));
fft_sig = 10*log10(abs(fft_sig));
fft_pu1 = 10*log10(abs(fft_pu1));
fft_pu2 = 10*log10(abs(fft_pu2));

df = 1/del_t;
f=linspace(0,df,numel(fft_sig));

fft_sig_axis = zeros(numel(f),1);
fft_pu1_axis = zeros(numel(f),1);
fft_pu2_axis = zeros(numel(f),1);
for ii=1:numel(f)
    fft_sig_axis(ii,1) = (c/(c/lam_sig + f(1,ii) - 8e11));
    fft_pu1_axis(ii,1) = (c/(c/lam_pu1 + f(1,ii) - 8e11));
    fft_pu2_axis(ii,1) = (c/(c/lam_pu2 + f(1,ii) - 8e11));
end

plot(fft_sig_axis,fft_sig)
hold on
plot(fft_pu1_axis,fft_pu1)
plot(fft_pu2_axis,fft_pu2)
plot(ASE_spectral_axis,ASE_spectral)
hold off

nonlinear_coefficient = gamma * div_n

toc
