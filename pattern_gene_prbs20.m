%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRBS(N=20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prbs = zeros(1,1048576);
regi = zeros(1,20);
regi(1,20) = 1;
f = zeros(1,20);
f(1,1)=1;
f(1,8)=1;

for i = 1:1048576
    prbs(1,i)=regi(1,1);
    buf = f(1,1)*regi(1,1);

    for k = 2:20
        buf = xor(buf,f(1,k)*regi(1,k));
    end

    for k = 2:20
        regi(1,k-1) = regi(1,k);
    end

    regi(1,20) = buf;
end

pattern = [zeros(1,5),prbs];


%%N=17
%{
prbs = zeros(1,131072);
regi = zeros(1,17);
regi(1,17) = 1;
f = zeros(1,17);
f(1,1)=1;
f(1,8)=1;

for i = 1:131072
    prbs(1,i)=regi(1,1);
    buf = f(1,1)*regi(1,1);

    for k = 2:17
        buf = xor(buf,f(1,k)*regi(1,k));
    end

    for k = 2:17
        regi(1,k-1) = regi(1,k);
    end

    regi(1,17) = buf;
end

pattern = [zeros(1,5),prbs];

%}

%pattern = [prbs,prbs];
%save-binary pattern_PRBS6 pattern