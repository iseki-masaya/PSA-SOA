%1/2*piˆÊ‘Š•ÏXŠÖ”
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_2 = phase_rev(E_1)

F_in = fft(E_1);
E_in_phase = angle(F_in);
E_in_power = abs(F_in);

E_in_phase2 = E_in_phase - pi/2;
E_2 = ifft(E_in_power .* exp(1i*E_in_phase2))/sqrt(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%