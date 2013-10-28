function [] = plot_tsignal( E_in,t )
%PLOT_TSIGNAL Summary of this function goes here
%   Detailed explanation goes here

subplot(2,1,1)
plot(t*1e9-0.1,abs(E_in).^2*10^3);
xlabel('time[ns]');
ylabel('Intensity[mW]');
xlim([0,0.5]);

subplot(2,1,2)
plot(t*1e9-0.1,unwrap(angle(E_in)));
xlabel('time[ns]');
ylabel('Phase[rad]');
xlim([0,0.5]);
ylim([-pi,pi])

end

