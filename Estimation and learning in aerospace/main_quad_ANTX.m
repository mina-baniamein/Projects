%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANT-X SIMULATOR - MAIN                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;
addpath('datasets','common','common/simulator-toolbox','common/simulator-toolbox/attitude_library','common/simulator-toolbox/trajectory_library');
clc;

%% Model parameters

% Initial model (state: longitudinal velocity, pitch rate, pitch angle; input: normalised pitching moment; outputs: state and longitudinal acceleration)

Xu=-0.1068;

Xq=0.1192;

Mu=-5.9755;

Mq=-2.6478;

Xd=-10.1647;

Md=450.71;

A=[Xu, Xq, -9.81; Mu, Mq, 0; 0, 1, 0];

B=[Xd; Md; 0];

C=[1, 0, 0; 0, 1, 0; 0, 0, 1; Xu, Xq, 0]; 

D=[0; 0; 0; Xd];

% Noise

%noise.Enabler = 0;
noise.Enabler = 1;

noise.pos_stand_dev = noise.Enabler * 0.0011;                            	%[m]

noise.vel_stand_dev = noise.Enabler * 0.01;                               %[m/s]

noise.attitude_stand_dev = noise.Enabler * deg2rad(0.33);                 %[rad]
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(1);                   %[rad/s]

% Delays

delay.position_filter = 1;
delay.attitude_filter = 1;
delay.mixer = 1;

%% Load controller parameters

parameters_controller                    

%% M injection example (sweeep: first column time vector, secondo column time history of pitching moment) 
% Selection on the input
% New
x = 0;
while (x ~= 1) && (x ~= 2)
    x = input('Select 1 for ExitationM and 2 for input u3211 : ');
    if x == 1
        load ExcitationM
    elseif x == 2
        Ts_3ord = 0.001; % Sampling time for u3211
        Amplitude = 0.1; % Amplitude
        Ripetition = 19;
        ExcitationM = u3211(Amplitude,Ts_3ord,Ripetition); % Function to generate u3211 input
    else
        disp('Any valid number was selected');
    end
end

% Plotting input signal
figure('Name','Input plot')
t_3ord=ExcitationM(:,1);
u_3ord=ExcitationM(:,2);
plot(t_3ord, u_3ord)
grid
title('Input Signal')
ylabel('u'); ylim([min(u_3ord) - 0.05, max(u_3ord) + 0.05]); xlabel('time [s]')

SetPoint=[0,0];

%% Values selected
t=ExcitationM(:,1);

simulation_time=t(end)-t(1);

%% Launch SIMULATOR

%sim Simulator_Single_Axis
simulation_data = sim('Simulator_Single_Axis','SrcWorkspace', 'current');

%% Cross-Validation Separation

u_3ord = simulation_data.Mtot; 
y_q = simulation_data.q;
y_ax = simulation_data.ax;
t_3ord=(0:sample_time:simulation_time)';
% Percentage
p_cross = 0.75;
%
u_cross = u_3ord(ceil(length(u_3ord)*p_cross)+1:end);
u_3ord = u_3ord(1:ceil(length(u_3ord)*p_cross));
y_q_cross = y_q(ceil(length(y_q)*p_cross)+1:end);
y_q = y_q(1:ceil(length(y_q)*p_cross));
y_ax_cross = y_ax(ceil(length(y_ax)*p_cross)+1:end);
y_ax = y_ax(1:ceil(length(y_ax)*p_cross));
t_3ord_cross = t_3ord(ceil(length(t_3ord)*p_cross)+1:end);
t_3ord = t_3ord(1:ceil(length(t_3ord)*p_cross));

figure
subplot(2,2,1)
hold on
plot(t_3ord, y_q,'r');
plot(t_3ord_cross, y_q_cross,'k'); grid on;
title('Input Signal')
ylabel('y_q'); ylim([min([min(y_q),min(y_q_cross)]) - 0.05, max([max(y_q),max(y_q_cross)]) + 0.05]); xlabel('time [s]')

subplot(2,2,2)
hold on
plot(t_3ord, y_ax,'r');
plot(t_3ord_cross, y_ax_cross,'k'); grid on;
title('Input Signal')
ylabel('y_a_x'); ylim([min([min(y_ax),min(y_ax_cross)]) - 0.05, max([max(y_ax),max(y_ax_cross)]) + 0.05]); xlabel('time [s]')

subplot(2,1,2)
hold on
plot(t_3ord, u_3ord,'r');
plot(t_3ord_cross, u_cross,'k'); grid on;
title('Input Signal')
ylabel('u'); ylim([min([min(u_3ord),min(u_cross)]) - 0.05, max([max(u_3ord),max(u_cross)]) + 0.05]); xlabel('time [s]')
%% CONTINUA CROSS VALIDATION
% Modifica bettering solution per metterlo prima e poi vediamo
%% TASK 1 : Grey-box model identification
% Model ID by grey-box model
% Input Assigned
y_grey = [y_q y_ax]; 
% Model call
model_fun = @drone_model_grey;
% Guess vector
guess = [Xu Xq Mu Mq Xd Md] + 0.1 * (-0.5 + rand(1,6)).*[Xu Xq Mu Mq Xd Md];
% Real Parameter Vector
real_parameters = [Xu; Xq; Mu; Mq; Xd; Md];
% Function call
[identification_grey, error_grey] = Model_identification(u_3ord,y_grey,guess,sample_time,model_fun,real_parameters);
A_grey = identification_grey.matrix{1}; B_grey = identification_grey.matrix{2};
C_grey = identification_grey.matrix{3}; D_grey = identification_grey.matrix{4};
%% TASK 2 : Black-box model identification
% Input Assigned
y_3ord = y_q;
Ts_3ord = sample_time;
% Spectral Figure and determination of the cutting frequency
spectral_analysis(y_3ord,Ts_3ord);
fc = input('Select the cutting frequency of the output : ');
[u_3ord,y_3ord] = bettering_solution(u_3ord,y_3ord,fc,Ts_3ord);
%% Definition of the better p value for system from output error
n = 3; % Order of the system
min_p = 2*n; max_p = 40; 
p = better_p(u_3ord,y_3ord,min_p,max_p,3,Ts_3ord,t_3ord);
disp('The minimum error on the output is on p :'); disp(p);
%% PBSID 
[D_PBSID,S_PBSID,V_PBSID,Y,N] = pbsid_1part(u_3ord,y_3ord(:,1),p);
n = input("looking the graph choose the order detected by the PBSID : ");
[A_PBSID,B_PBSID,C_PBSID,K_PBSID] = pbsid_2part(D_PBSID,n,S_PBSID,V_PBSID,Y,N,p,u_3ord,y_3ord);
%% Structuring with greyest
y_id = lsim(ss(A_PBSID,B_PBSID,C_PBSID,D_PBSID), u_3ord, t_3ord);
theta0 = [Xu Xq Mu Mq Xd Md] + 0.01 * (-0.5 + rand(1,6)).*[Xu Xq Mu Mq Xd Md];
model_fun = @drone_model_PBSID;
[identification_PBSID, error_PBSID] = greyest_structuring(u_3ord,y_id,Ts_3ord,theta0,model_fun,real_parameters);
% Graphic comparison of outputs
figure ("Name","Graphic comparison of real outputs wrt identified sys")
plot(t_3ord,y_3ord,'k',t_3ord,lsim(ss(A_PBSID,B_PBSID,C_PBSID,D_PBSID,Ts_3ord), u_3ord, t_3ord),'r');
[A_PBSID,B_PBSID,C_PBSID,D_PBSID] = drone_model_grey(identification_PBSID.parameters,0);
%% Analysis on freqresp
w = logspace(-3, 3, 1000);
%
sys_real = ss(A,B,C([2,4],:),D([2,4],:));
sys_PBSID = ss(A_PBSID,B_PBSID,C_PBSID,D_PBSID);
sys_grey = ss(A_grey,B_grey,C_grey,D_grey);
G_real = tf(sys_real);
G_PBSID = tf(sys_PBSID);
G_grey = tf(sys_grey);
% Calcolo della risposta in frequenza
H_real = freqresp(sys_real, w);
H_PBSID = freqresp(sys_PBSID, w);
H_grey = freqresp(sys_grey, w);
% H è un array 3D: per sistemi SISO, lo "squeeze" lo rende vettore 1D
H_real = squeeze(H_real);
H_PBSID = squeeze(H_PBSID);
H_grey = squeeze(H_grey);

% Modulo e fase
magnitude_real = abs(H_real);
phase_real = angle(H_real);
magnitude_PBSID = abs(H_PBSID);
phase_PBSID = angle(H_PBSID);
magnitude_grey = abs(H_grey);
phase_grey = angle(H_grey);

% Plot
figure('Name','Errori');
hold on
subplot(2,2,1);
semilogx(w, (20*log10(magnitude_real)-20*log10(magnitude_PBSID)));grid on;
title('Modulo (dB)');
xlabel('Frequenza (rad/s)');
ylabel('|H(jω)| [dB]');

subplot(2,2,2);
semilogx(w, (20*log10(magnitude_real)-20*log10(magnitude_grey))); grid on;
title('Modulo (dB)');
xlabel('Frequenza (rad/s)');
ylabel('|H(jω)| [dB]');

subplot(2,2,3);
semilogx(w, (rad2deg(phase_real)-rad2deg(phase_PBSID)));grid on;
title('Fase (gradi)');
xlabel('Frequenza (rad/s)');
ylabel('Fase [°]');

subplot(2,2,4);
semilogx(w,(rad2deg(phase_real)-rad2deg(phase_grey))); grid on;
title('Fase (gradi)');
xlabel('Frequenza (rad/s)');
ylabel('Fase [°]');

error_PBSID
error_grey
%% Cross-Validation
% PBSID
A = A_PBSID; B = B_PBSID;
C([2,4],:) = C_PBSID; D([2,4],:) = D_PBSID;
simulation_data_PBSID = sim('Simulator_Single_Axis','SrcWorkspace', 'current');
%%
%plot((0:sample_time:simulation_time)',simulation_data_PBSID.q,'k',(0:sample_time:simulation_time)',simulation_data.q,'r--')
clear x
err_PBSID_Cross = simulation_data_PBSID.q - simulation_data.q;
max_err_Cross = max(abs(err_PBSID_Cross));

figure('Name', 'Final Validation on Unseen Data');
subplot(2, 1, 1);
hold on;
plot((0:sample_time:simulation_time)', simulation_data.q, 'b-', 'LineWidth', 1.5);
plot((0:sample_time:simulation_time)', simulation_data_PBSID.q, 'r--', 'LineWidth', 1.5);
hold off;
grid on;
title('Model Performance on Validation Set','FontSize',15);
xlabel('Time [s]');
ylabel('Pitch Rate q [rad/s]');
legend('True Output', 'Predicted Output');
subplot(2, 1, 2);
plot((0:sample_time:simulation_time)', err_PBSID_Cross, 'k-');
grid on;
title('Prediction Error (Predicted - True)','FontSize',15);
xlabel('Time [s]');
ylabel('Error [rad/s]');
ylim([-max_err_Cross*1.2, max_err_Cross*1.2]); 
%% Delete temporary files

if exist('slprj','dir')
    rmdir('slprj', 's')                                                    
end

%% END OF CODE

