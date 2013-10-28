num = size(A_sig,1);

for i = 1:1:num-div_n
    E_out_sig(i,1)=A_sig(i+div_n,div_n);
end

Generate_Constellation(E_out_sig,bit_d);
Generate_Eyepattern(abs(E_out_sig).^2,bit_d);

[output_sig,Q]=DPSK_receiver(E_out_sig,bit_d);
sigma = get_sigma(E_out_sig,bit_d)
Q
