clear; clc; close all;

%% 1 - Define a simple first-order discrete-time SISO system
Ts = 0.1;   % Sampling time [s]

A_true = 0.8;    % scalar
B_true = 0.5;    % scalar
C_true = 1;      % scalar
D_true = 0;      % scalar

sys_true = ss(A_true,B_true,C_true,D_true,Ts);

%% 2 - Input: Chirp excitation
Tfinal = 50;
t = (0:Ts:Tfinal-Ts)';
u = chirp(t, 0.05, Tfinal, 1.0);   % Frequency sweep
y = lsim(sys_true,u,t);

%% 3 - PBSID pipeline (simplified for test case)
n = 1;   % true order
min_p = 2*n; max_p = 20; p = 10;  % prediction horizon

[D,S,V,Ys,N] = pbsid_1part(u,y,p);
[A_est,B_est,C_est,K_est] = pbsid_2part(D,n,S,V,Ys,N,p,u,y);

sys_pbsid = ss(A_est,B_est,C_est,D,Ts);
y_pbsid = lsim(sys_pbsid,u,t);

%% 4 - Grey-box structuring with idgrey + greyest
% Pack true parameters as a vector
real_parameters = [A_true; B_true; C_true; D_true];

% Initial guess from PBSID + noise
theta0 = [A_est; B_est; C_est; D] + 0.01*randn(numel(real_parameters),1);

% Create grey-box model definition (discrete-time, sample time = Ts)
sys_init = idgrey(@siso_model_PBSID, theta0, 'd', Ts);

% Identification using time-domain iddata
sim_data = iddata(y,u,Ts);
options = greyestOptions('Display','on','SearchMethod','lsqnonlin');
options.SearchOptions.FunctionTolerance = 1e-6;
estimated_model = greyest(sim_data, sys_init, options);

% Extract results
theta_id = estimated_model.Report.Parameters.ParVector;
fit_id   = estimated_model.Report.Fit.FitPercent;

[A_grey,B_grey,C_grey,D_grey] = siso_model_PBSID(theta_id,Ts);

sys_grey = ss(A_grey,B_grey,C_grey,D_grey,Ts);
y_grey = lsim(sys_grey,u,t);

%% 5 - Plots
figure;
plot(t,y,'k','LineWidth',1.5); hold on
plot(t,y_pbsid,'r--','LineWidth',1.2)
plot(t,y_grey,'b-.','LineWidth',1.2)
legend('True','PBSID','Grey-Structured')
xlabel('Time [s]'); ylabel('Output'); grid on

%% 6 - Report
disp('--- True system ---');  
disp(A_true); 
disp(B_true); 
disp(C_true); 
disp(D_true); 
disp('B*C ='); 
disp(B_true*C_true);

disp('--- PBSID estimate ---'); 
disp(A_est); 
disp(B_est); 
disp(C_est); 
disp(D); 
disp('B*C ='); 
disp(B_est*C_est);

disp('--- Grey-structured estimate ---'); 
disp(A_grey); 
disp(B_grey); 
disp(C_grey); 
disp(D_grey); 
disp('B*C ='); 
disp(B_grey*C_grey);

disp(['Greyest Fit [%]: ', num2str(fit_id)]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Helper model function for idgrey (first order SISO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B,C,D] = siso_model_PBSID(theta, Ts, ~)
    % theta contains [a; b; c; d]
    A = theta(1);
    B = theta(2);
    C = theta(3);
    D = theta(4);
end


