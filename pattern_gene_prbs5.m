%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRBS(N=5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prbs = zeros(1,31);
regi = zeros(1,5);
regi(1,5) = 1;
f = zeros(1,5);
f(1,1)=1;
f(1,2)=1;

for i = 1:31
    prbs(1,i)=regi(1,1);
    buf = f(1,1)*regi(1,1);

    for k = 2:5
        buf = xor(buf,f(1,k)*regi(1,k));
    end

    for k = 2:5
        regi(1,k-1) = regi(1,k);
    end

    regi(1,5) = buf;
end

pattern = [zeros(1,5),prbs,zeros(1,5)];
%pattern = [prbs,prbs];

%save-binary pattern_PRBS6 pattern