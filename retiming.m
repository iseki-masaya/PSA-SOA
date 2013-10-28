function [ E_out ] = retiming( E_in,bit_d )

num = size(E_in,1);
point_div = 1;
E_out = [];

for i = 1:1:num/bit_d
    for ii = -point_div:1:point_div
        E_out = [E_out;E_in(i * bit_d -bit_d/2 + ii,1)];
    end    
end

%   Detailed explanation goes here


end

