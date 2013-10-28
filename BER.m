

term;                       %室温
eta_i           = 0.8;      %外部量子効率
n_sp            = 4;        %自然放出係数
m_t             = 2;        %モード数
lambda          = 1500*10^(-9); % 波長
B_e;                        % 信号帯域
dat_len;
Res             = 50;        %外部抵抗[Ω]
B_opt           = 100*10^9;  %BPFの帯域
L_0             = -7;
L_1             = -16;
L_2             = -3;
L_3             = -13;
G_0             = 7;
G_1             = 16;
G_2             = 6;
NF              = 2;
G               = 20;
V_th = 0;
P_in = input_level_sig;
SIG = (L_0*G_0)*(L_1*G_1)*(L_2*G_2)*L_3;
ASE = (G_0-1)*(L_1*G_1)*(L_2*G_2)*L_3+(G_1-1)*(L_2*G_2)*L_3+(G_2-1)*L_3;

sigma_sshot = 2 * e * (e* eta_i*B*P_in/(h*ω)*SIG;
sigma_spshot = 2*e*(e*eta_i*B*n_sp*m_t*B_opt)*ASE;
sigma_ssp = 4*(e*eta_i*P_in/(h*ω)*(e*eta_i*n_sp*B)*SIG*ASE;
sigma_spsp = 4*(e*eta_i*n_sp)^2*m_t*B_opt*B*ASE^2;
sigma_th;
sigma_1 =sigma_sshot+sigma_spshot+sigma_ssp+sigma_spsp+sigma_th;
sigma_0 =+sigma_spshot+sigma_spsp+sigma_th;
sigma = (sigma_1+sigma_0)/2;
 func1 = @(x) x+x.^2;
>> quad(func1,0,1)
ans =
    0.8333
>> quadl(func1,0,1)
ans =
    0.8333
func1 = @(x) 1/(sqrt(Pi)*sigma)*(1/dat_len * exp(-(x-e*eta_i*SIG*G*Res*P_in/(h*ω))/(sqrt(2)*sigma*Res/sqrt(2));
ans1 = quad(func1,-100,V_th);    
func2 = @(x) (dat_len-1)/dat_len * exp(-(x/(sqrt(2)*sigma*Res/sqrt(2))^2);
ans2 = quad(func2,V_th,100);
BER = ans1 + ans2;








