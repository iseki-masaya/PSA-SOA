
function [complex_waveform, Power, Phase] = Fiber_Transmission(input_waveform,dt, input_phase, D, S, L, wavelength, center_wavelength)
parameter;
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

% -------------------------------------------
n2 = 1.22*10^(-14);	
Step = 100;		%?½?½Ô•ï¿½?½?½?½?½
dL = L/Step;
dL2 = dL/2;	

%Dw = -(wavelength)^2*D*L/2.99792458e8;
Dw = -(wavelength)^2*D*dL2/2.99792458e8;

%Sw = (wavelength^3/(2.99792458e8)^2)*(2*D*L+wavelength*S*L);
Sw = (wavelength^3/(2.99792458e8)^2)*(2*D*dL2+wavelength*S*dL2);

E_in = sqrt(input_waveform) .* exp(1i*input_phase);

%F_in = fftshift(fft(E_in));
F_in = fftshift(fft(E_in));

fs = 1/dt;
% dt=del_t
df = fs/dat_size_P;

n = 0:1:(dat_size_P-1);

%f = (n*df - fs/2)';
f = (n*df- fs/2)';
%?½?½?½S?½?½g?½?½?½?½0?½?½


% t = 0:del_t:del_t*(dat_len-1);
% dat_size_P=dat_len
f0 = 2.99792458e8/center_wavelength;

fc = 2.99792458e8/wavelength;

% H_fiber = exp(-1i * (0.5*Dw*(f+fc-f0).^2 + (1/6)*Sw*(f+fc-f0).^3));
% H_fiber = exp(-1i *(0.5*Dw*f.^2 + (1/6)*Sw*(f.^3 + 3*f.^2*(fc-f0))));
% H_fiber = exp(-1i *(0.5*Dw*f.^2));
H_fiber = exp(-1i * (0.5*Dw*f.^2 + (1/6)*Sw*(f.^3 + 3*f.^2*(fc-f0))));
%plot(f,angle(H_fiber));
% xlabel('Frequency (Hz)')
% ylabel('Phase(rad)')
% saveas(gcf,'F_phase_spectol.jpg');
% xlim([0,3*10^10]);

for s = 0:1:Step-1
    F_out = F_in .* H_fiber;
    E_out = ifft(ifftshift(F_out));
    
    E_out = signal_amp(E_out,-0.2*dL/1e3);
    
    E_out = sf(E_out,wavelength,n2,dL);
    F_in2 = fftshift(fft(E_out));
    F_out2 = F_in2 .* H_fiber;
    F_in = F_out2;
    

end


F_out = F_in;
E_out = ifft(ifftshift(F_out));
E_out = signal_amp(E_out,20);

complex_waveform = E_out;
% Phase    = unwrap(angle(E_out) + pi*0.5*((-1).^(n+1)+1)');
Power = (abs(E_out)).^2;

Phase    = (angle(E_out));

%{
F_out = F_in .* H_fiber;

E_out = ifft(ifftshift(F_out));
complex_waveform = E_out;
% Phase    = unwrap(angle(E_out) + pi*0.5*((-1).^(n+1)+1)');
Power = (abs(E_out)).^2;

Phase    = (angle(E_out));
%}

%Phase    = unwrap(angle(E_out) + pi*0.5*((-1).^(n+1)+1)');

% Phase    = (angle(E_out));


