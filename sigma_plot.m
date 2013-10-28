%4pi
%sigma_phase = [pi()/20,pi()/18,pi()/16,pi()/14,pi()/12,pi()/10,pi()/8,pi()/6,pi()/4];
%sigma_phai_in =[2.3080,2.4796,2.5904,2.6649,2.7180,2.7575,2.7880,2.8123,2.8320];
sigma_phai_in =[2.8320,2.8123,2.7880,2.7575,2.7180,2.6649,2.5904,2.4796,2.3080];
%sigma_phai_out = [2.2764,2.5573,2.5999,2.6301,2.6554,2.6731,2.6831,2.6901,2.6958];
sigma_phai_out = [2.6958,2.6901,2.6831,2.6731,2.6554,2.6301,2.5999,2.5573,2.2764];

% plot(sigma_phai_in,sigma_phai_out);
% xlabel('Input sigma');
% ylabel('Output sigma');
% xlim([2.3,3]);
% ylim([2.3,3]);
grid
grid on
plot(sigma_phai_out./sigma_phai_in);
xlabel('Phase[rad]');
ylabel('ratio of input sigma of phase and Output sigma of phase');
%set(gca,'XTick',[pi()/20,pi()/18,pi()/16,pi()/14,pi()/12,pi()/10,pi()/8,pi()/6,pi()/4])
set(gca,'XTickLabel',{'pi()/20','pi()/18','pi()/16','pi()/14','pi()/12','pi()/10','pi()/8','pi()/6','pi()/4'})
% xlim([2.3,3]);