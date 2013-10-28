function [ out_phase ] = angledef( E_in )
E_in_phase = angle(E_in);
% roop = 1;
% 
% while roop
%     
% if max(E_in_phase)>0
%     E_in_phase = E_in_phase-2*pi();
%     if (max(E_in_phase))<0
%         roop = 0;
%     end
% end
% 
% if max(E_in_phase)<0
%     E_in_phase=E_in_phase + 2*pi();
%     if (max(E_in_phase))>0
%         roop = 0;
%     end
% end
% end
out_phase = E_in_phase;
% 
% %ANGLEDEF Summary of this function goes here
% %   Detailed explanation goes here


end

