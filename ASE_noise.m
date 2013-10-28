function [ P_a ] = ASE_noise( n )
%ASE_NOISE Summary of this function goes here
%   Detailed explanation goes here
%mu = 20;
%mu = 30;
mu = 1;
h               = 6.6260693e-34; 
c               = 2.99792458e8;
lambda          = 1550*10^(-9); 

B = 10 * 10^9;
G = (100-1)*n;
P_a = mu * (h * c/lambda)*(G)*B;

end
