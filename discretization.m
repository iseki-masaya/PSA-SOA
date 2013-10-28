function [ H_func_out ] = discretization(H_func_in,Lcos_res,Lcos_grad,signal_bandwidth )

%LCOS�̉𑜓x�ɂ�鐧��
parameter;
H_size = size(H_func_in,1);
df = 1/del_t;
H_size_eff = round(H_size*signal_bandwidth/df);

H_func_in(H_size_eff:1:H_size-H_size_eff) = 0;

Lcos_delta = round(H_size/Lcos_res*df/signal_bandwidth/2);  

for index = 1:Lcos_delta:H_size
        for roop = 1:1:Lcos_delta-1
                H_func_in(index+roop,1)=H_func_in(index,1);
        end
end
%�T�C�Y����
if size(H_func_in,1)-H_size~= 0
        H_func_in=H_func_in(1:H_size,1);
end

Lcos_power_tmp=abs(H_func_in);
Lcos_phase_tmp=angle(H_func_in);

%LCOS�̊K�����ɂ�鐧��


for i =1:1:H_size
    
        for ii = 1:1:Lcos_grad+1          
                if Lcos_phase_tmp(i,1)>0
                    if Lcos_phase_tmp(i,1) < 2*pi()/Lcos_grad*ii
                        Lcos_phase_tmp(i,1) = 2*pi()/Lcos_grad*(ii-1);
                        break
                    end
                else
                    if Lcos_phase_tmp(i,1) > -2*pi()/Lcos_grad*ii
                        Lcos_phase_tmp(i,1) = -2*pi()/Lcos_grad*(ii-1);
                        break
                    end
                end
        end
end          
H_func_out = Lcos_power_tmp .* exp(1i*Lcos_phase_tmp);

%�T�C�Y����
if size(H_func_in,1)-H_size ~= 0
        H_func_in=H_func_in(1:H_size,1);
end



end

