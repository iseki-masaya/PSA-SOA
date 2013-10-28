function [ I_sat_out ] = I_sat( d )
eta_i = 0.8;      %�ʎq����
A = 3.0 * 10e-15;       %�z���f�ʐ�[cm^2]
lambda          = 1550*10^(-9); % �g��
h               = 6.6260693e-34;            %?v?????N????[J/s]
c               = 2.99792458e8;             %����[m/s]
alpha_0 = 1500;         %�s�O�a�z����[cm^-1] 
tau = 10*10e-12;
I_sat_out = (eta_i * A * tau/(h * c/lambda * alpha_0 * d))^(-1);

%2cm���炢


end

