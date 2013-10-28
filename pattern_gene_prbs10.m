%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRBS(N=9)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prbs = zeros(1,1023);
regi = zeros(1,10);
regi(1,10) = 1;
f = zeros(1,10);
f(1,1)=1;
f(1,8)=1;

for i = 1:1023
    prbs(1,i)=regi(1,1);
    buf = f(1,1)*regi(1,1);

    for k = 2:10
        buf = xor(buf,f(1,k)*regi(1,k));
    end

    for k = 2:10
        regi(1,k-1) = regi(1,k);
    end

    regi(1,10) = buf;
end

pattern = [zeros(1,5),prbs];
%pattern = [prbs,prbs];

%save-binary pattern_PRBS6 pattern