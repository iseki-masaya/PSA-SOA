[E_out, Power_out, Phase_out] = SPM(P_input_sig,del_t,Phase_input_sig,D_SSMF,S_SSMF,100e3,1550e-9, 1565e-9); 
                                                                                                                                  
%subplot(2,1,1); plot(t,E_MZI_bar_out)
%subplot(2,1,2); plot(t,E_MZI_cross_out)