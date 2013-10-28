function [complex_waveform, Power, Phase] = SPM(input_waveform,dt, input_phase, D, S, L, wavelength, center_wavelength)
          

% data size verify
% ------------------------------------------
dat_size_P_r = size(input_waveform,1);
dat_size_P_c = size(input_waveform,2);

if dat_size_P_r == 1
    if dat_size_P_c == 1
        error('signal data is invalid !!');
    end
    dat_size_P = dat_size_P_c;
    input_waveform = input_waveform';
elseif dat_size_P_c == 1
    dat_size_P = dat_size_P_r;
end

dat_size_phase_r = size(input_phase,1);
dat_size_phase_c = size(input_phase,2);

if dat_size_phase_r == 1
    if dat_size_phase_c == 1
        error('signal data is invalid !!');
    end
    dat_size_phase = dat_size_phase_c;
    input_phase = input_phase';
elseif dat_size_phase_c == 1
    dat_size_phase = dat_size_phase_r;
end

if dat_size_P ~= dat_size_phase
    error('signal data is invalid !!');
end

t = get_time(dat_size_P);

% -------------------------------------------

speed           = 10e9;                     %?????M?????r?b?g???[?g
bit_d           = 160;                       %1?r?b?g?f?[?^????????
del_t           = 1 / (speed * bit_d);   
n2 = 1.22*10^(-14);		%�K���X���`�萔[m/V]%????????[s]      
Step = 100;		%��ԕ�����
dL = L/Step;
dL2 = dL/2;		%�v�Z��̋�ԕ�����

f0 = 2.99792458e8/center_wavelength;
fc = 2.99792458e8/wavelength;

Dw = -(wavelength)^2*(D)*dL2/2.99792458e8;
Sw = (wavelength^3/(2.99792458e8)^2)*(2*D*dL2+wavelength*S*dL2);

E_in = sqrt(input_waveform) .* exp(1i*input_phase);
fs = 1/dt;
df = fs/dat_size_P;
n = 0:1:(dat_size_P-1);
f = (n*df- fs/2)';
f0 = 2.99792458e8/center_wavelength;

fc = 2.99792458e8/wavelength;

H_fiber = exp(-1i *(0.5*Dw*f.^2)+(1/6)*Sw*(f.^3 + 3*f.^2*(fc-f0)));
% H_fiber = exp(-1i *(0.5*Dw*f.^2));

E_out = E_in;
for s = 0:1:Step-1
	F_in = fftshift(fft(E_in));
	F_out = F_in .* H_fiber;
	E_out = ifft(ifftshift(F_out));

    E_out = sf(E_out,wavelength,n2,dL);
	F_in = fftshift(fft(E_out));
	F_out = F_in .* H_fiber;
	E_in =  ifft(ifftshift(F_out));
%     E_in = signal_amp(E_in,-0.2*dL/1e3);

%         subplot(2,1,2)
%          plot(real(E_in));
end
% E_in = signal_amp(E_in,20);
E_out = E_in;

t1 = t * 10^9;
% E_in_dat = E_in .* sqrt(10^3);
% E_in2_dat = E_in2 .* sqrt(10^3);
% hold on
% plot(t1,abs(E_in_dat).^2,'r')
% xlabel('time[ns]');
% ylabel('Intensity[mW]');
% plot(t1,abs(E_in2_dat).^2)
% hold off
complex_waveform = E_out;
Power = (abs(E_out)).^2;
Phase    = (angle(E_out));
% Phase    = unwrap(angle(E_out) + pi*0.5*((-1).^(n+1)+1)');
%Phase    = unwrap(angle(E_out) + pi*0.5*((-1).^(n+1)+1)');

% Phase    = (angle(E_out));
end

