A_I_data=csvread('waveform6.csv');
A_Q_data=csvread('waveform7.csv');
A_I=A_I_data(17:4000,2);
A_Q=A_Q_data(17:4000,2);
A=(A_I+i*A_Q);
Generate_Constellation(A,40);

E=0.01*exp(i*0);
D(1,1)=A(1,1)-E;

for C=2:1:size(A_I)/20
    D(C,1)=A((C-1)*20,1)-D(C-1,1);
end

F=abs(D);
G=angle(D);
