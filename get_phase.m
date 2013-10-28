function [E_phase] = get_phase(E_in,P_input_sig)

E_base = sqrt(P_input_sig);
E_phase = angle(E_in./E_base);