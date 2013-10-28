function [ana] = Generate_Eyepattern(E_in,bit_d);


 E_in = E_in/max(abs(E_in));

 t = get_time(size(E_in,1));
 
% F_in = fft(E_in);
% i = size(F_in,2)
% for ii = 1:1:i
%     if abs(F_in(1,ii).^2) < 20000;
%         F_in(1,ii) = 0;
%      end
% end
% plot(abs(F_in).^2);
% E_in = ifft(F_in);


%Fs =10*10^9*bit_d ;
Fs =40*10^9*bit_d ;
Rs = size(E_in,1)/bit_d;
nSamps = bit_d/2;
%nSamps = bit_d;
rollOff = 0.5;

hMod = modem.pskmod('M', 2, 'PhaseOffset', 0);

hFig = figure(40);
%plot(t,E_in);
%managescattereyefig(hFig);
%{
eyeObj = commscope.eyediagram(...
    'SamplingFrequency', Fs, ...
    'MinimumAmplitude', -1, ...
    'MaximumAmplitude', 1, ...
    'SamplesPerSymbol', nSamps, ...
    'PlotType', '2D Color', ... 
    'ColorScale', 'linear',...
    'OperationMode', 'real signal')
%}
eyeObj = commscope.eyediagram(...
    'SamplingFrequency', Fs, ...
    'MinimumAmplitude', -1, ...
    'MaximumAmplitude', 1, ...
    'NumberOfStoredTraces' , 1000, ... 
    'SamplesPerSymbol', nSamps);
%eyeObj. PlotType = '2D Line';

%analyze(eyeObj)
%eyeObj.PlotTimeOffset = 1/10^(9)*80;
%eyeObj.update(0.5*E_in);
eyeObj.update(E_in);
managescattereyefig(hFig, eyeObj, 'right');
cmap = jet(64);
cmap(1,:) = [0 0 0];
plot(eyeObj, cmap);
%plot(eyeObj , jet(64));
%plot(eyeObj);
%eyeObj
save 'eyediagram' ,exportdata(eyeObj);

%plot(eyeObj, cmap);
% analyze(eyeObj)


