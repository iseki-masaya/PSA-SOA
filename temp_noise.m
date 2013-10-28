rmatrix = randn(10000000,1);    %100ƒR
% mean1 = mean(rmatrix)
% disp(     '=0')
% var1 = var(rmatrix)
% disp('    =1')

num_mean = 0
num_var = 4

rmatrix2 = rmatrix * sqrt(num_var) + num_mean;
rmatrix2 = rmatrix2/2;
mean2 = mean(rmatrix2)
disp('      =1.5751e-006')
var2 = var(rmatrix2)
('    1')
hist(rmatrix);

mean3 = mean(rmatrix2.^2)
var3 = var(rmatrix2.^2)
disp('---------------------------------------------------------')