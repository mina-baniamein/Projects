function [p] = better_p(u,y,min_p,max_p,n,Ts,t,step)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation and Learning in Aerospace Project A.Y. 24-25 
% Function determining the better value of p in order to get less erros on
% the output funtion

% Outputs : p : [1 x 1] p value with less disturbances

% Authors:  Alessandro Castelanelli Oddo (alessandro.castelanelli@polimi.it)
%            (@polimi.it)                     
%            (@polimi.it)                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Array of p and inizialization of errors array
if nargin<8
    p_values = min_p:1:max_p; 
else
    p_values = min_p:step:max_p; 
end
error = zeros(size(p_values));
% Linear Search on PBSID
for i = 1:length(p_values)
    [D_id, S_svd, V_svd, Y, N] = pbsid_1part(u, y(:,1), p_values(i),1);
    [A_id, B_id, C_id, ~] = pbsid_2part(D_id, n, S_svd, V_svd, Y, N, p_values(i), u, y);
    sys_id = ss(A_id, B_id, C_id, D_id, Ts);
    error(i) = norm(y - lsim(sys_id, u, t)); % Errore di simulazione
end
% Plotting errors
figure('Name','y_output error by p')
semilogy(p_values, error, 'o-')
xlabel('p')
ylabel('Simulation Error')
title('P evolution error','FontSize',20)
grid on
% Inset on t = 10 h
axes('Position', [0.6, 0.6, 0.3, 0.3]);
hold on; grid on; box on;
title('Inset around p = [33, 36]')
semilogy(p_values, error, 'o-')
xlim([33, 36]);

% P search for the minimum error
p = p_values(find(error==min(error)));
clc
disp(min(error))
end