%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRBS generation (N=6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prbs = zeros(1,63);
regi = zeros(1,6);
regi(1,6) = 1;
f = zeros(1,6);
f(1,1)=1;
f(1,2)=1;

for i = 1:63
    prbs(1,i)=regi(1,1);
    buf = f(1,1)*regi(1,1);

    for k = 2:6
        buf = xor(buf,f(1,k)*regi(1,k));
    end

    for k = 2:6
        regi(1,k-1) = regi(1,k);
    end

    regi(1,6) = buf;
end

pattern = [prbs,prbs];
