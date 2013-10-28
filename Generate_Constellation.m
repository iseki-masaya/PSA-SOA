function [binary] = Generate_Constellation(E_in,Bit_d);


hMod = modem.pskmod('M', 2, 'PhaseOffset', 0);
Rup = Bit_d;  % up sampling rate
sym = E_in;
cap_E = [];
for i_cap = 1:1:size(E_in)/Bit_d
    cap_E = [cap_E;abs(E_in(Bit_d * i_cap - Bit_d/2,1))];
end
sym = sym./mean(cap_E);       %?½?½?½K?½?½
%Fs =10;
speed=40e9;
Fs = speed*Bit_d;

sym = circshift(sym,[Rup/2,1]);

%figure(60)
plot(real(sym));
hScope = commscope.ScatterPlot(...
    'SamplingFrequency' , Fs,...
    'SamplesPerSymbol' , Rup);
hScope.SymbolRate;
hScope.Constellation = hMod.Constellation;

update(hScope, sym)
hScope.PlotSettings.Constellation = 'on';
hScope.PlotSettings.SignalTrajectory = 'on';
autoscale(hScope);

binary = 1;
end