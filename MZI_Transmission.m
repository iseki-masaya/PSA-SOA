%��͂͂Ƃ肠������͐M��
% close all
tic
% main;
Q_map=[];
sigma_map=[];
% load('make_colormap0.mat','Q_map','sigma_map');
Gaussnoise = 0;
if_fluc=0;
if_sa = 3;
if_sa_i =2;

% for  i_alpha = 0.02:0.02:2.4
%     for i_i_s = 0.01*10^(-4):0.01*10^(-4):1.2*10^(-4);
i_alpha = 0.04;
i_i_s =0.2*10^(-4);
% i_alpha = 1.2;
% i_i_s =0.6*10^(-4);


% for d = 0.37:0.01:4
% for P_i = -30:0.5:30
[E_MZI_bar_out, Power_MZI_out, Phase_MZI_out,Q_map,sigma_map] = MZI_Filter(Q_map,...
    sigma_map,i_alpha,i_i_s,...
    P_input_sig,del_t,Phase_input_sig,D_SSMF,S_SSMF,100e3,1550e-9,...
    1550e-9,Gaussnoise,if_sa,if_fluc,if_sa_i); 
%      end
%      save('make_colormap0.mat','Q_map','sigma_map');
%      toc
%  end

    %�K�E�X�G��
%subplot(2,1,1); plot(t,E_MZI_bar_out)
%subplot(2,1,2); plot(t,E_MZI_cross_out)