function [u3211_vect] = u3211(A,Ts,rip)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation and Learning in Aerospace Project A.Y. 24-25 
% Function creating a  long sequence u3211 with amplitude as "A" at sample time
% "Ts" in a total time of 8 with an ammount of N samples 

% Outputs : u3211 : [N x 2] Sequence u3211 [time,input]

% Authors:  Alessandro Castelanelli Oddo (alessandro.castelanelli@polimi.it)
%            (@polimi.it)                     
%            (@polimi.it)                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_1 = ones(1/Ts,1);
t_2 = ones(2/Ts,1);
t_3 = ones(3/Ts,1);

u_section = [0*t_2;A*t_3;-A*t_2;A*t_1;-A*t_1];
n = length(u_section);
u = [0*t_3;u_section];
for i = 1:rip
    u = [u; u_section];
end
u = [u;0*t_3];

t = (0:Ts:length(u)*Ts)'; t = t(1:(end-1),1);
% % Plot using stairs
% figure;
% plot(t, u, 'b');
% xlabel('Time [s]');
% ylabel('u');
% grid on;

u3211_vect = [t,u];

end




