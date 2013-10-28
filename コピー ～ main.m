%clear                               %initialize

disp('     ===================================     ')
disp('     =====    Dual-pump SOA PSA    =====     ')
disp('     =====            DPSK         =====     ')
disp('     =====        ver. 1.00        =====     ')
disp('     ===================================     ')
                                    %display the version of this program
                                    
%load 'input_sig'

parameter                           %setting of parameter (constants)

input_condition                     %setting of input conditions

waveform_generation          %making optical waveform

tic                                 %start time count

TMM                                 %Transfer Matrix Method

toc                                 %stop time count

num=size(input_sig,1);

for i=1:1:num-42
    P_for_sig_shift(i,1)=P_for_sig(i+41,div_n-1);
end

%{
output_phase = zeros(73,1);
output_power = zeros(73,1);

for r=1:73
   
    output_phase(r,1) = angle(P_for_sig((780+(40*r)),41)) - angle(P_for_sig(390,41));
    output_power(r,1) = abs(P_for_sig((780+(40*r)),41));
    
end
%}
%{
%output
phi=unwrap(phi_for_sig(:,div_n)-phi_for_sig(:,1));
P_out=P_for_sig(:,div_n-1);

%phase matching condition
phase_condition=2*pi*(2/lam_sig - 1/lam_pu1 - 1/lam_pu2)+...
    gamma*div_n*10^(-3)*(10^(input_level_pu1/10)+10^(input_level_pu2/10)-10^(input_level_sig/10));
%}
%CSV file export
%{
csvwrite('P_sig_out.csv',P_out);
csvwrite('phi_sig_out.csv',phi(1:dat_len-1,:));
csvwrite('P_pu1_out.csv',P_for_pu1(:,div_n-1));
csvwrite('P_pu2_out.csv',P_for_pu2(:,div_n-1));

csvwrite('P_spm_at_41th.csv',P_spm_sig(:,div_n-1));
csvwrite('P_xpm_at_41th.csv',P_xpm_sig(:,div_n-1));
csvwrite('P_fwm_at_41th.csv',P_fwm_sig(:,div_n-1));
csvwrite('phi_spm_at_41th.csv',phi_spm_sig(:,div_n-1));
csvwrite('phi_xpm_at_41th.csv',phi_xpm_sig(:,div_n-1));
csvwrite('phi_fwm_at_41th.csv',phi_fwm_sig(:,div_n-1));
csvwrite('P_sig_at_41th_without_nonlinear.csv',abs((F_sig(:,41)).^2));
csvwrite('phi_sig_at_41th_without_nonlinear.csv',angle(F_sig(:,41)));
csvwrite('phi_relationship_between_sig_and_spm_at_41th.csv',unwrap(angle(F_sig(:,41))-angle(spm_sig(:,41))));
csvwrite('phi_relationship_between_sig_and_xpm_at_41th.csv',unwrap(angle(F_sig(:,41))-angle(xpm_sig(:,41))));
csvwrite('phi_relationship_between_sig_and_fwm_at_41th.csv',unwrap(angle(F_sig(:,41))-angle(fwm_sig(:,41))));

csvwrite('P_spm_at_2th.csv',P_spm_sig(:,2));
csvwrite('P_xpm_at_2th.csv',P_xpm_sig(:,2));
csvwrite('P_fwm_at_2th.csv',P_fwm_sig(:,2));
csvwrite('P_sig_at_2th_without_nonlinear.csv',abs((F_sig(:,2)).^2));
csvwrite('phi_sig_at_2th_without_nonlinear.csv',angle(F_sig(:,2)));
csvwrite('phi_relationship_between_sig_and_spm_at_2th.csv',unwrap(angle(F_sig(:,2))-angle(spm_sig(:,2))));
csvwrite('phi_relationship_between_sig_and_xpm_at_2th.csv',unwrap(angle(F_sig(:,2))-angle(xpm_sig(:,2))));
csvwrite('phi_relationship_between_sig_and_fwm_at_2th.csv',unwrap(angle(F_sig(:,2))-angle(fwm_sig(:,2))));
%}