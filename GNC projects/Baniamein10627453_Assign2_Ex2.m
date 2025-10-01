%Spacecraft Guidance and Navigation
%Assignment #2
%Author: Mina Baniamein

%% Point 1
clear ; close all; clc;
%% Kernel unload
cspice_kclear(); %unload previous kernels section
%% Meta kernel
cspice_furnsh('kernels\assignment2.txt'); 
%%
% Constant for arcseconds to radians conversions
arcsec2rad = pi / (180*3600);
typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)
% TLE SMOS
longstr1 = '1 36036U 09059A   24323.76060260  .00000600  00000-0  20543-3 0  9995';
longstr2 = '2 36036  98.4396 148.4689 0001262  95.1025 265.0307 14.39727995790658';
% Initialize the satrec structure, using the function twoline2rv
satrec = twoline2rv(longstr1, longstr2, typerun,'e', opsmode, whichconst);
% Get TLE epoch
[year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str);

fprintf('Satellite num ID: %d\n', satrec.satnum);
fprintf('TLE reference epoch: UTC %s', sat_epoch_str); %t ref= UTC 2024-11-18 18:15:16.064640 
% Evaluate the TLE
[satrec,rteme,vteme] = sgp4(satrec,  0.0);
% Get the osculating orbital elements 
elts = cspice_oscelt( [rteme;vteme], sat_epoch_et, satrec.mu );

fprintf('*** Osculating orbital elements ***')
fprintf('SMA   [km]:  %.5f', elts(1)/(1-elts(2)));
fprintf('ECC   [km]:  %.8f', elts(2));
fprintf('INC  [deg]: %.5f', elts(3)*cspice_dpr()); %conversion from rad to degrees
fprintf('RAAN [deg]: %.5f', elts(4)*cspice_dpr());
fprintf('ARGP [deg]: %.5f', elts(5)*cspice_dpr());
fprintf('M.AN [deg]: %.5f', elts(6)*cspice_dpr());

a = elts(1)/(1-elts(2));

ddpsi = -0.115178;
ddeps = -0.007395;

et0 = cspice_str2et(sat_epoch_str); % Epoch we want to perform our conversion TEME to ECI, expressed in seconds past J2000 epoch
ttt = cspice_unitim(et0, 'ET', 'TDT')/cspice_jyear()/100; % %number of centuries past J2000 for TEMEtoECI conversion 

% Transform TEME to ECI vectors
ateme = [0;0;0];
[reci0, veci0, aeci0] = teme2eci(rteme, vteme, ateme, ttt, ddpsi, ddeps); %reci0 and veci0 are 3x1 column
fprintf(1,' reci0 %14.7f %14.7f %14.7f',reci0 );
fprintf(1,' veci0 %14.9f %14.9f %14.9f',veci0 );
fprintf(1,' aeci0 %14.9f %14.9f %14.9f',aeci0 );

mu = satrec.mu;

%Keplerian propagation of ECI states
ti = 0;
tf = 2*pi*sqrt((a^3)/mu); %orbital period, a in norma è uguale ad aeci
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[~, xx] = ode113(@(t,x) keplerian_rhs(t, x, mu), [ti tf], [reci0; veci0], options);

sat_epoch0_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [2024,11,18,20,30,0]);
sat_epochf_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [2024,11,18,22,15,0]);
et_i=cspice_str2et(sat_epoch0_str); %expressed in seconds past J2000 epoch
et_f=cspice_str2et(sat_epochf_str); %expressed in seconds past J2000 epoch

% Equally spaced vector, with one minute time-step FOR KOUROU STATION AND
% SVALBARD
npoints = round((et_f-et_i)/60.0)+1;    %minutes past J2000
et_vec = linspace(et_i, et_f, npoints); %minutes sized spacing

%Equally spaced vector, with 30s time-step because sampling time of
%TROLL
npoints_troll = round((et_f-et_i)/30.0)+1; %minutes past J2000
et_vec_troll = linspace(et_i, et_f, npoints_troll); %half minutes sized spacing 

%Computation of ECI state at 18/11/24 20:30
[~, xxi] = ode113(@(t,x) keplerian_rhs(t, x, mu), [0 (et_i-et0)], [reci0; veci0], options);
reci_i=xxi(end,1:3)';
veci_i=xxi(end,4:6)';
%Computation of ECI states from at 18/11/24 to 20:30 to 22:15
[~, xxf] = ode113(@(t,x) keplerian_rhs(t, x, mu), et_vec, [reci_i; veci_i], options);
[~, xxf_troll] = ode113(@(t,x) keplerian_rhs(t, x, mu), et_vec_troll, [reci_i; veci_i], options);

reci_f=xxf(end,1:3)';
reci_f_troll=xxf_troll(end,1:3)'; %check that the final state computed is the same with a refined time grid
veci_f=xxf(end,4:6)';  %probabilmente inutile questo stato rv alla etf


% Compute antenna angles, satellite range and range-rate of every station
[SMOS_azimuth_kourou, SMOS_elevation_kourou, SMOS_range_kourou, SMOS_range_rate_kourou] = ...
    antenna_pointing('KOUROU', et_vec, [xxf(:,1:3)';xxf(:,4:6)']);

[SMOS_azimuth_troll, SMOS_elevation_troll, SMOS_range_troll, SMOS_range_rate_troll] = ...
    antenna_pointing('TROLL', et_vec_troll, [xxf_troll(:,1:3)';xxf_troll(:,4:6)']);

[SMOS_azimuth_svalbard, SMOS_elevation_svalbard, SMOS_range_svalbard, SMOS_range_rate_svalbard] = ...
    antenna_pointing('SVALBARD', et_vec, [xxf(:,1:3)';xxf(:,4:6)']);

%azimuth plot
figure(2)
plot(et_vec/cspice_spd(), SMOS_azimuth_kourou*cspice_dpr(),'b')
hold on
plot(et_vec_troll/cspice_spd(), SMOS_azimuth_troll*cspice_dpr(),'r')
plot(et_vec/cspice_spd(), SMOS_azimuth_svalbard*cspice_dpr(),'g')
title('SMOS Azimuth')
xlabel('Epoch [MJD2000]')
ylabel('Azimuth [deg]')
legend('KOUROU','TROLL','SVALBARD')
%elevation plot
figure(3)
plot(et_vec/cspice_spd(), SMOS_elevation_kourou*cspice_dpr(),'b')
hold on
plot(et_vec_troll/cspice_spd(), SMOS_elevation_troll*cspice_dpr(),'r')
plot(et_vec/cspice_spd(), SMOS_elevation_svalbard*cspice_dpr(),'g')
title('SMOS Elevation')
xlabel('Epoch [MJD2000]')
ylabel('Elevation [deg]')
legend('KOUROU','TROLL','SVALBARD')

% Plot passes of all stations (azimuth and elevation)
% Elevation
figure(4)
i_visibility_kourou = SMOS_elevation_kourou > deg2rad(6);
plot(et_vec(i_visibility_kourou)./cspice_spd(), SMOS_elevation_kourou(i_visibility_kourou).*cspice_dpr(),'*')
hold on
grid on
i_visibility_troll = SMOS_elevation_troll > deg2rad(0);
plot(et_vec_troll(i_visibility_troll)./cspice_spd(), SMOS_elevation_troll(i_visibility_troll).*cspice_dpr(),'*')
i_visibility_svalbard = SMOS_elevation_svalbard > deg2rad(8);
plot(et_vec(i_visibility_svalbard)./cspice_spd(), SMOS_elevation_svalbard(i_visibility_svalbard).*cspice_dpr(),'*')
axis([et_vec(1)/cspice_spd(),et_vec(end)/cspice_spd(),0, 90])
xlabel('Epoch [MJD2000]')
ylabel('Elevation [deg]')
legend('KOUROU','TROLL','SVALBARD')
title('Elevation during visibiity windows')

% Azimuth
figure
i_visibility_kourou = SMOS_elevation_kourou > deg2rad(6);
plot(et_vec(i_visibility_kourou)./cspice_spd(), SMOS_azimuth_kourou(i_visibility_kourou).*cspice_dpr(),'*')
hold on
grid on
i_visibility_troll = SMOS_elevation_troll > deg2rad(0);
plot(et_vec_troll(i_visibility_troll)./cspice_spd(), SMOS_azimuth_troll(i_visibility_troll).*cspice_dpr(),'*')
i_visibility_svalbard = SMOS_elevation_svalbard > deg2rad(8);
plot(et_vec(i_visibility_svalbard)./cspice_spd(), SMOS_azimuth_svalbard(i_visibility_svalbard).*cspice_dpr(),'*')
% axis([et_vec(1)/cspice_spd(),et_vec(end)/cspice_spd(),0, 90])
xlabel('Epoch [MJD2000]')
ylabel('Azimuth [deg]')
legend('KOUROU','TROLL','SVALBARD')
title('Azimuth during visibiity windows')

clc

TW_kourou = et_vec(i_visibility_kourou);
Kourou_TW_start = cspice_et2utc(TW_kourou(1) ,'C', 3);
Kourou_TW_ending = cspice_et2utc(TW_kourou(end) ,'C', 3);
disp(['Kourou visibiity window beginning: ', num2str(Kourou_TW_start)]);
disp(['Kourou visibiity window ending: ', num2str(Kourou_TW_ending)]);
disp(' ')



TW_troll = et_vec_troll(i_visibility_troll);
Troll_TW_start = cspice_et2utc(TW_troll(1) ,'C', 3);
Troll_TW_ending = cspice_et2utc(TW_troll(end) ,'C', 3);
disp(['Troll visibiity window beginning: ', num2str(Troll_TW_start)]);
disp(['Troll visibiity window ending: ', num2str(Troll_TW_ending)]);
disp(' ')


TW_svalbard = et_vec(i_visibility_svalbard);
Svalbard_TW_start = cspice_et2utc(TW_svalbard(1) ,'C', 3);
Svalbard_TW_ending = cspice_et2utc(TW_svalbard(end) ,'C', 3);
disp(['Svalbard visibiity window beginning: ', num2str(Svalbard_TW_start)]);
disp(['Svalbard visibiity window ending: ', num2str(Svalbard_TW_ending)]);
disp(' ')


%% Point 2 
% Kourou measuraments computation
stationName = 'KOUROU';
[range_real_kourou, azimuth_real_kourou, elevation_real_kourou] = measuraments(TW_kourou,et0,stationName,satrec,ddpsi,ddeps);

figure
hold on
grid on
i_visibility_kourou_meas = elevation_real_kourou > deg2rad(6);
plot(TW_kourou(i_visibility_kourou_meas)./cspice_spd(), elevation_real_kourou(i_visibility_kourou_meas).*cspice_dpr(),'*')
plot(et_vec(i_visibility_kourou)./cspice_spd(), SMOS_elevation_kourou(i_visibility_kourou).*cspice_dpr(),'*')
legend('Measuared Elevation','Predicted Elevation')
title('Kourou measuraments')

TW_kourou = TW_kourou(i_visibility_kourou_meas);
Kourou_TW_start = cspice_et2utc(TW_kourou(1) ,'C', 3);
Kourou_TW_ending = cspice_et2utc(TW_kourou(end) ,'C', 3);

disp(['Measured Kourou visibiity window beginning: ', num2str(Kourou_TW_start)]);
disp(['Measured Kourou visibiity window ending: ', num2str(Kourou_TW_ending)]);
disp(' ')

% Troll measuraments computation
stationName = 'TROLL';
[range_real_troll, azimuth_real_troll, elevation_real_troll] = measuraments(TW_troll,et0,stationName,satrec,ddpsi,ddeps);

figure
hold on
grid on
i_visibility_troll_meas = elevation_real_troll > deg2rad(0);
plot(TW_troll(i_visibility_troll_meas)./cspice_spd(), elevation_real_troll(i_visibility_troll_meas).*cspice_dpr(),'*')
plot(et_vec_troll(i_visibility_troll)./cspice_spd(), SMOS_elevation_troll(i_visibility_troll).*cspice_dpr(),'*')
legend('Measuared Elevation','Predicted Elevation')
title('Troll measuraments')

TW_troll = TW_troll(i_visibility_troll_meas);
Troll_TW_start = cspice_et2utc(TW_troll(1) ,'C', 3);
Troll_TW_ending = cspice_et2utc(TW_troll(end) ,'C', 3);

disp(['Measured Troll visibiity window beginning: ', num2str(Troll_TW_start)]);
disp(['Measured Troll visibiity window ending: ', num2str(Troll_TW_ending)]);
disp(' ')

% Svalbard measuraments computation
stationName = 'SVALBARD';
[range_real_svalbard, azimuth_real_svalbard, elevation_real_svalbard] = measuraments(TW_svalbard,et0,stationName,satrec,ddpsi,ddeps);

figure
hold on
grid on
i_visibility_svalbard_meas = elevation_real_svalbard > deg2rad(8);
plot(TW_svalbard(i_visibility_svalbard_meas)./cspice_spd(), elevation_real_svalbard(i_visibility_svalbard_meas).*cspice_dpr(),'*')
plot(et_vec(i_visibility_svalbard)./cspice_spd(), SMOS_elevation_svalbard(i_visibility_svalbard).*cspice_dpr(),'*')
legend('Measuared Elevation','Predicted Elevation')
title('Svalbard measuraments')

TW_svalbard = TW_svalbard(i_visibility_svalbard_meas);
Svalbard_TW_start = cspice_et2utc(TW_svalbard(1) ,'C', 3);
Svalbard_TW_ending = cspice_et2utc(TW_svalbard(end) ,'C', 3);

disp(['Svalbard visibiity window beginning: ', num2str(Svalbard_TW_start)]);
disp(['Svalbard visibiity window ending: ', num2str(Svalbard_TW_ending)]);
disp(' ')

%% Point 3
%Real measuraments storing for lsqnonlin
sigma_range = 0.01;                %[km]
sigma_azimuth = deg2rad(0.125);    %[rad]
sigma_elevation = deg2rad(0.125);  %[rad]
%Weigths matrix
meas_noise_cov = diag([sigma_range^2;sigma_azimuth^2;sigma_elevation^2]);
W_m = inv(sqrtm(meas_noise_cov));

%% a) Kourou only, pure Keplerian motion
clc
%Matrixes inizialization for storing
meas_real = zeros(length(range_real_kourou),3);
meas_real(:,1) = range_real_kourou;
meas_real(:,2) = azimuth_real_kourou;
meas_real(:,3) = elevation_real_kourou;

meas_real = {meas_real};
TW =  {TW_kourou};
stationNames = {'KOUROU'};
% Keplerian dynamics choice
fundynamics = @(t,x) keplerian_rhs(t,x,mu);
% Cost function
fun = @(x) costfunction(x, TW, W_m, meas_real, stationNames, fundynamics, et_i); %input: initial guess, 6x1 vector

%Guess is predicted ECI state  at 20:30 with keplerian propagation
x0_a = [reci_i(:,1); veci_i(:,1)]; 


% Call lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
options.FunctionTolerance = 1e-10;
options.StepTolerance = 1e-10;

[x_a,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun, x0_a , [], [], options);
%%
% Print results
disp('x - x_0 =');
disp(x_a - x0_a);
norm(x_a(1:3) - x0_a(1:3))
norm(x_a(4:6) - x0_a(4:6))

% Covariance computation
Jac_a = full(jacobian);
P_ls_a = resnorm/(length(residual)-length(x_a)).*inv(Jac_a.'*Jac_a);
trace_r_a = sqrt(trace(P_ls_a(1:3,1:3)));
trace_v_a = sqrt(trace(P_ls_a(4:6,4:6)));

% Keplerian covariance computation
J_kep_a = computeJacobianKeplerian(x_a, mu);
P_kep_a = J_kep_a*P_ls_a*J_kep_a';

%% b) All stations, pure Keplerian motion
clc
meas_real_kourou = zeros(length(range_real_kourou),3);
meas_real_kourou(:,1) = range_real_kourou;
meas_real_kourou(:,2) = azimuth_real_kourou;
meas_real_kourou(:,3) = elevation_real_kourou;

meas_real_troll = zeros(length(range_real_troll),3);
meas_real_troll(:,1) = range_real_troll;
meas_real_troll(:,2) = azimuth_real_troll;
meas_real_troll(:,3) = elevation_real_troll;

meas_real_svalbard = zeros(length(range_real_svalbard),3);
meas_real_svalbard(:,1) = range_real_svalbard;
meas_real_svalbard(:,2) = azimuth_real_svalbard;
meas_real_svalbard(:,3) = elevation_real_svalbard;

meas_real = {meas_real_kourou, meas_real_troll, meas_real_svalbard};

stationNames = {'KOUROU','TROLL','SVALBARD'};
TW = {TW_kourou,TW_troll,TW_svalbard};

% Keplerian dynamics choice
fundynamics = @(t,x) keplerian_rhs(t,x,mu);

% Cost function
fun = @(x) costfunction(x, TW, W_m, meas_real, stationNames, fundynamics, et_i); %input: initial guess, 6x1 vector

%Guess is predicted ECI state  at 20:30 with keplerian propagation
x0_a = [reci_i(:,1); veci_i(:,1)]; 

% Call lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
options.FunctionTolerance = 1e-10;
options.StepTolerance = 1e-10;

[x_b,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun, x0_a , [], [], options);
%%
% Print results
disp('x - x_0 =');
disp(x_b - x0_a);
norm(x_b(1:3) - x0_a(1:3))
norm(x_b(4:6) - x0_a(4:6))

% Covariance computation
Jac_b = full(jacobian);
P_ls_b = resnorm/(length(residual)-length(x_b)).*inv(Jac_b.'*Jac_b);
trace_r_b = sqrt(trace(P_ls_b(1:3,1:3)));
trace_v_b = sqrt(trace(P_ls_b(4:6,4:6)));

% Keplerian covariance computation
J_kep_b = computeJacobianKeplerian(x_b, mu);
P_kep_b = J_kep_b*P_ls_b*J_kep_b';

%% c) All stations, J2 perturbed keplerian motion
clc
meas_real_kourou = zeros(length(range_real_kourou),3);
meas_real_kourou(:,1) = range_real_kourou;
meas_real_kourou(:,2) = azimuth_real_kourou;
meas_real_kourou(:,3) = elevation_real_kourou;

meas_real_troll = zeros(length(range_real_troll),3);
meas_real_troll(:,1) = range_real_troll;
meas_real_troll(:,2) = azimuth_real_troll;
meas_real_troll(:,3) = elevation_real_troll;

meas_real_svalbard = zeros(length(range_real_svalbard),3);
meas_real_svalbard(:,1) = range_real_svalbard;
meas_real_svalbard(:,2) = azimuth_real_svalbard;
meas_real_svalbard(:,3) = elevation_real_svalbard;

meas_real = {meas_real_kourou, meas_real_troll, meas_real_svalbard};

stationNames = {'KOUROU','TROLL','SVALBARD'};
TW = {TW_kourou,TW_troll,TW_svalbard};

% Keplerian dynamics choice
fundynamics = @(t,x) J2_rhs(t,x,mu);

% Cost function
fun = @(x) costfunction(x, TW, W_m, meas_real, stationNames, fundynamics, et_i); %input: initial guess, 6x1 vector

%Guess is predicted ECI state  at 20:30 with keplerian propagation
x0_a = [reci_i(:,1); veci_i(:,1)]; 

% Call lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
options.FunctionTolerance = 1e-10;
options.StepTolerance = 1e-10;

[x_c,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun, x0_a , [], [], options);
%%
% Print results
disp('x - x_0 =');
disp(x_c - x0_a);
norm(x_c(1:3) - x0_a(1:3))
norm(x_c(4:6) - x0_a(4:6))

% Covariance computation
Jac_c = full(jacobian);
P_ls_c = resnorm/(length(residual)-length(x_c)).*inv(Jac_c.'*Jac_c);
trace_r_c = sqrt(trace(P_ls_c(1:3,1:3)));
trace_v_c = sqrt(trace(P_ls_c(4:6,4:6)));

% Keplerian covariance computation
J_kep_c = computeJacobianKeplerian(x_c, mu);
P_kep_c = J_kep_c*P_ls_c*J_kep_c';


%% Point 4
% KOUROU - TROLL
clc
meas_real_kourou = zeros(length(range_real_kourou),3);
meas_real_kourou(:,1) = range_real_kourou;
meas_real_kourou(:,2) = azimuth_real_kourou;
meas_real_kourou(:,3) = elevation_real_kourou;

meas_real_troll = zeros(length(range_real_troll),3);
meas_real_troll(:,1) = range_real_troll;
meas_real_troll(:,2) = azimuth_real_troll;
meas_real_troll(:,3) = elevation_real_troll;

meas_real = {meas_real_kourou, meas_real_troll};

stationNames = {'KOUROU','TROLL'};
TW = {TW_kourou,TW_troll};

% Keplerian dynamics choice
fundynamics = @(t,x) J2_rhs(t,x,mu);

% Cost function
fun = @(x) costfunction(x, TW, W_m, meas_real, stationNames, fundynamics, et_i); %input: initial guess, 6x1 vector

%Guess is predicted ECI state  at 20:30 with keplerian propagation
x0_a = [reci_i(:,1); veci_i(:,1)]; 

% Call lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
options.FunctionTolerance = 1e-10;
options.StepTolerance = 1e-10;

[x_KT,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun, x0_a , [], [], options);

% Print results
disp('x - x_0 =');
disp(x_KT - x0_a);

% Covariance computation
Jac_KT = full(jacobian);
P_ls_KT = resnorm/(length(residual)-length(x_KT)).*inv(Jac_KT.'*Jac_KT);
trace_r_KT = sqrt(trace(P_ls_KT(1:3,1:3)));
trace_v_KT = sqrt(trace(P_ls_KT(4:6,4:6)));
disp(['R trace = ', num2str(trace_r_KT), '; V trace = ', num2str(trace_v_KT)])

% Keplerian covariance computation
J_kep_KT = computeJacobianKeplerian(x_KT, mu);
P_kep_KT = J_kep_KT*P_ls_KT*J_kep_KT';
%%
clc
sigma_a_i_KT = [P_kep_KT(1,1),P_kep_KT(3,3)];
disp(['Deviation on a = ', num2str(sigma_a_i_KT(1)), '; Deviation on i = ', num2str(sigma_a_i_KT(2))])

%% KOUROU - SVALBARD (best)

meas_real_kourou = zeros(length(range_real_kourou),3);
meas_real_kourou(:,1) = range_real_kourou;
meas_real_kourou(:,2) = azimuth_real_kourou;
meas_real_kourou(:,3) = elevation_real_kourou;

meas_real_svalbard = zeros(length(range_real_svalbard),3);
meas_real_svalbard(:,1) = range_real_svalbard;
meas_real_svalbard(:,2) = azimuth_real_svalbard;
meas_real_svalbard(:,3) = elevation_real_svalbard;

meas_real = {meas_real_kourou, meas_real_svalbard};

stationNames = {'KOUROU','SVALBARD'};
TW = {TW_kourou,TW_svalbard};

% Keplerian dynamics choice
fundynamics = @(t,x) J2_rhs(t,x,mu);

% Cost function
fun = @(x) costfunction(x, TW, W_m, meas_real, stationNames, fundynamics, et_i); %input: initial guess, 6x1 vector

%Guess is predicted ECI state  at 20:30 with keplerian propagation
x0_a = [reci_i(:,1); veci_i(:,1)]; 

% Call lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
options.FunctionTolerance = 1e-10;
options.StepTolerance = 1e-10;

[x_KS,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun, x0_a , [], [], options);


% Covariance computation
Jac_KS = full(jacobian);
P_ls_KS = resnorm/(length(residual)-length(x_KS)).*inv(Jac_KS.'*Jac_KS);
trace_r_KS = sqrt(trace(P_ls_KS(1:3,1:3)));
trace_v_KS = sqrt(trace(P_ls_KS(4:6,4:6)));
disp(['R trace = ', num2str(trace_r_KS), '; V trace = ', num2str(trace_v_KS)])

% Keplerian covariance computation
J_kep_KS = computeJacobianKeplerian(x_KS, mu);
P_kep_KS = J_kep_KS*P_ls_KS*J_kep_KS';


sigma_a_i_KS = [P_kep_KS(1,1),P_kep_KS(3,3)];
disp(['Deviation on a = ', num2str(sigma_a_i_KS(1)), '; Deviation on i = ', num2str(sigma_a_i_KS(2))])

%% TROLL - SVALBARD

meas_real_troll = zeros(length(range_real_troll),3);
meas_real_troll(:,1) = range_real_troll;
meas_real_troll(:,2) = azimuth_real_troll;
meas_real_troll(:,3) = elevation_real_troll;

meas_real_svalbard = zeros(length(range_real_svalbard),3);
meas_real_svalbard(:,1) = range_real_svalbard;
meas_real_svalbard(:,2) = azimuth_real_svalbard;
meas_real_svalbard(:,3) = elevation_real_svalbard;

meas_real = {meas_real_troll, meas_real_svalbard};

stationNames = {'TROLL','SVALBARD'};
TW = {TW_troll,TW_svalbard};

% Keplerian dynamics choice
fundynamics = @(t,x) J2_rhs(t,x,mu);

% Cost function
fun = @(x) costfunction(x, TW, W_m, meas_real, stationNames, fundynamics, et_i); %input: initial guess, 6x1 vector

%Guess is predicted ECI state  at 20:30 with keplerian propagation
x0_a = [reci_i(:,1); veci_i(:,1)]; 

% Call lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
options.FunctionTolerance = 1e-10;
options.StepTolerance = 1e-10;

[x_ST,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun, x0_a , [], [], options);

% Print results
disp('x - x_0 =');
disp(x_ST - x0_a);

% Covariance computation
Jac_ST = full(jacobian);
P_ls_ST = resnorm/(length(residual)-length(x_ST)).*inv(Jac_ST.'*Jac_ST);
trace_r_ST = sqrt(trace(P_ls_ST(1:3,1:3)));
trace_v_ST = sqrt(trace(P_ls_ST(4:6,4:6)));
disp(['R trace = ', num2str(trace_r_ST), '; V trace = ', num2str(trace_v_ST)])

% Keplerian covariance computation
J_kep_ST = computeJacobianKeplerian(x_ST, mu);
P_kep_ST = J_kep_ST*P_ls_ST*J_kep_ST';


sigma_a_i_ST = [P_kep_ST(1,1),P_kep_ST(3,3)];
disp(['Deviation on a = ', num2str(sigma_a_i_ST(1)), '; Deviation on i = ', num2str(sigma_a_i_ST(2))])

%%
% KOUROU - TROLL

clc
meas_real_kourou = zeros(length(range_real_kourou),3);
meas_real_kourou(:,1) = range_real_kourou;
meas_real_kourou(:,2) = azimuth_real_kourou;
meas_real_kourou(:,3) = elevation_real_kourou;

meas_real_troll = zeros(length(range_real_troll),3);
meas_real_troll(:,1) = range_real_troll;
meas_real_troll(:,2) = azimuth_real_troll;
meas_real_troll(:,3) = elevation_real_troll;

meas_real = {meas_real_kourou, meas_real_troll};

stationNames = {'KOUROU','TROLL'};
TW = {TW_kourou,TW_troll};

% Keplerian dynamics choice
fundynamics = @(t,x) J2_rhs(t,x,mu);

% Cost function
fun = @(x) costfunction(x, TW, W_m, meas_real, stationNames, fundynamics, et_i); %input: initial guess, 6x1 vector

%Guess is predicted ECI state  at 20:30 with keplerian propagation
x0_a = [reci_i(:,1); veci_i(:,1)]; 

% Call lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
options.FunctionTolerance = 1e-10;
options.StepTolerance = 1e-10;

[x_KT,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun, x0_a , [], [], options);

% Print results
disp('x - x_0 =');
disp(x_KT - x0_a);

% Covariance computation
Jac_KT = full(jacobian);
P_ls_KT = resnorm/(length(residual)-length(x_KT)).*inv(Jac_KT.'*Jac_KT);
trace_r_KT = sqrt(trace(P_ls_KT(1:3,1:3)));
trace_v_KT = sqrt(trace(P_ls_KT(4:6,4:6)));

% Keplerian covariance computation
J_kep_KT = computeJacobianKeplerian(x_KT, mu);
P_kep_KT = J_kep_KT*P_ls_KT*J_kep_KT';

sigma_a_i_KT = [P_kep_KT(1,1),P_kep_KT(3,3)];
%% KOUROU - SVALBARD 
clc
meas_real_kourou = zeros(length(range_real_kourou),3);
meas_real_kourou(:,1) = range_real_kourou;
meas_real_kourou(:,2) = azimuth_real_kourou;
meas_real_kourou(:,3) = elevation_real_kourou;

meas_real_svalbard = zeros(length(range_real_svalbard),3);
meas_real_svalbard(:,1) = range_real_svalbard;
meas_real_svalbard(:,2) = azimuth_real_svalbard;
meas_real_svalbard(:,3) = elevation_real_svalbard;

meas_real = {meas_real_kourou, meas_real_svalbard};

stationNames = {'KOUROU','SVALBARD'};
TW = {TW_kourou,TW_svalbard};

% Keplerian dynamics choice
fundynamics = @(t,x) J2_rhs(t,x,mu);

% Cost function
fun = @(x) costfunction(x, TW, W_m, meas_real, stationNames, fundynamics, et_i); %input: initial guess, 6x1 vector

%Guess is predicted ECI state  at 20:30 with keplerian propagation
x0_a = [reci_i(:,1); veci_i(:,1)]; 

% Call lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
options.FunctionTolerance = 1e-10;
options.StepTolerance = 1e-10;

[x_KS,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun, x0_a , [], [], options);

% Print results
disp('x - x_0 =');
disp(x_KS - x0_a);

% Covariance computation
Jac_KS = full(jacobian);
P_ls_KS = resnorm/(length(residual)-length(x_KS)).*inv(Jac_KS.'*Jac_KS);
trace_r_KS = sqrt(trace(P_ls_KS(1:3,1:3)));
trace_v_KS = sqrt(trace(P_ls_KS(4:6,4:6)));

% Keplerian covariance computation
J_kep_KS = computeJacobianKeplerian(x_KS, mu);
P_kep_KS = J_kep_KS*P_ls_KS*J_kep_KS';

sigma_a_i_KS = [P_kep_KS(1,1),P_kep_KS(3,3)];

%% TROLL - SVALBARD (migliore per 1 giorno)
clc
meas_real_troll = zeros(length(range_real_troll),3);
meas_real_troll(:,1) = range_real_troll;
meas_real_troll(:,2) = azimuth_real_troll;
meas_real_troll(:,3) = elevation_real_troll;

meas_real_svalbard = zeros(length(range_real_svalbard),3);
meas_real_svalbard(:,1) = range_real_svalbard;
meas_real_svalbard(:,2) = azimuth_real_svalbard;
meas_real_svalbard(:,3) = elevation_real_svalbard;

meas_real = {meas_real_troll, meas_real_svalbard};

stationNames = {'TROLL','SVALBARD'};
TW = {TW_troll,TW_svalbard};

% Keplerian dynamics choice
fundynamics = @(t,x) J2_rhs(t,x,mu);

% Cost function
fun = @(x) costfunction(x, TW, W_m, meas_real, stationNames, fundynamics, et_i); %input: initial guess, 6x1 vector

%Guess is predicted ECI state  at 20:30 with keplerian propagation
x0_a = [reci_i(:,1); veci_i(:,1)]; 

% Call lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
options.FunctionTolerance = 1e-10;
options.StepTolerance = 1e-10;

[x_ST,resnorm,residual,exitflag,~,~,jacobian] = lsqnonlin(fun, x0_a , [], [], options);

% Print results
disp('x - x_0 =');
disp(x_ST - x0_a);

% Covariance computation
Jac_ST = full(jacobian);
P_ls_ST = resnorm/(length(residual)-length(x_ST)).*inv(Jac_ST.'*Jac_ST);
trace_r_ST = sqrt(trace(P_ls_ST(1:3,1:3)));
trace_v_ST = sqrt(trace(P_ls_ST(4:6,4:6)));

% Keplerian covariance computation
J_kep_ST = computeJacobianKeplerian(x_ST, mu);
P_kep_ST = J_kep_ST*P_ls_ST*J_kep_ST';

sigma_a_i_ST = [P_kep_ST(1,1),P_kep_ST(3,3)];





%% Function
function J = computeJacobianKeplerian(state, mu)
% DESCRIPTION
%   This function computes the numerical Jacobian matrix of the mapping 
%   from Cartesian state variables to Keplerian orbital elements using 
%   finite differences. Each state component is individually perturbed, 
%   and the resulting variations in the Keplerian elements are used to 
%   estimate the partial derivatives.
%
% INPUT
%   state   [6x1]   Cartesian state vector: [x; y; z; vx; vy; vz]                         [km, km/s]
%   mu      [1x1]   Gravitational parameter of the central body                          [km^3/s^2]
%
% OUTPUT
%   J       [6x6]   Jacobian matrix of Keplerian elements w.r.t. Cartesian state         [varies]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB
    
J = zeros(6,6); % Initialize Jacobian

% Perturbations
drx = perturbation(state(1));
dry = perturbation(state(2));
drz = perturbation(state(3));
dvx = perturbation(state(4));
dvy = perturbation(state(5));
dvz = perturbation(state(6));
pert = diag([drx, dry, drz, dvx, dvy, dvz]);

%Reference mean state kep elements
[a, e, i, Omega, omega, theta] = car2kep(state(1:3), state(4:6), mu);

for j=1:6
    state_p = state + pert(:,j);
    [a_p, e_p, i_p, Omega_p, omega_p, theta_p] = car2kep(state_p(1:3), state_p(4:6), mu);
    J(:,j) = ([a_p, e_p, i_p, Omega_p, omega_p, theta_p] - [a, e, i, Omega, omega, theta])./pert(j,j);
end

end

function dx = perturbation(x)
% DESCRIPTION
%   Generates a floating-point-based perturbation vector for Jacobian computation, 
%   scaled relative to the input state vector.
%
% INPUT
%   x    [Nx1]     State vector to perturb                                              [-]
%
% OUTPUT
%   dx   [Nx1]     Perturbation vector for each state component                         [-]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

dx = sqrt(eps)*max(abs(x),1);
end

function [Az,El,range,range_rate] = antenna_pointing(Station_Name,tt,xx)
% DESCRIPTION
%   Computes azimuth, elevation, range, and range rate for a ground station 
%   tracking a spacecraft, based on their relative positions and velocities.
%
% INPUT
%   Station_Name   [char]     Name of the ground station as defined in SPICE kernels         [-]
%   tt             [1xN]      Vector of ephemeris times                                      [s]
%   xx             [6xN]      Spacecraft state vector in ECI frame (position and velocity)   [km, km/s]
%
% OUTPUT
%   Az             [1xN]      Azimuth angles of the spacecraft relative to the station       [rad]
%   El             [1xN]      Elevation angles of the spacecraft relative to the station     [rad]
%   range          [1xN]      Range (distance) from station to spacecraft                    [km]
%   range_rate     [1xN]      Range rate (relative radial velocity)                          [km/s]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

    % Initialize vectors
    l = length(tt);
    Az = zeros(1,l); 
    El = zeros(1,l);
    range = zeros(1,l); 
    range_rate = zeros(1,l); 

    for i = 1 : l

        % Compute relative position in ECI
        xx_station_ECI = cspice_spkezr(Station_Name, tt(i), 'J2000', 'NONE', 'EARTH');
        xx_relative_ECI = xx(:,i) - xx_station_ECI;

        % Coversion from ECI to topocentric frame
        ROT_ECI2TOPO = cspice_sxform('J2000', [Station_Name, '_TOPO'], tt(i));
        xx_relative_topo = ROT_ECI2TOPO*xx_relative_ECI;

        % Compute range, range rate and azimuth and elevation
        [range(i), Az(i), El(i)] = cspice_reclat(xx_relative_topo(1:3));
        range_rate(i) = dot(xx_relative_topo(4:6),xx_relative_topo(1:3)/range(i));

    end

end

function residual_tot = costfunction(x0, TW, W_m, meas_real, stationNames, fundynamics, et_i) 
% DESCRIPTION
%   This function computes the total residual vector between real and predicted 
%   measurements (range, azimuth, elevation) for a satellite observed from multiple 
%   ground stations over specified time windows. It propagates the satellite state 
%   using provided dynamics, predicts measurements via antenna pointing, and computes 
%   weighted residuals using a measurement covariance weighting matrix.
%
% INPUT
%   x0            [6x1]       Initial guess for the satellite state vector                          [varies]
%   TW            {1xm}       Cell array containing time windows for each station                   [s]
%   W_m           [3x3]       Measurement noise covariance weighting matrix                         [-]
%   meas_real     {1xm}       Cell array containing real (simulated/noisy) measurements             [km, rad]
%   stationNames  {1xm}       Cell array containing ground station names as registered in SPICE     [-]
%   fundynamics   [function]  Handle to the RHS dynamics function used for state propagation        [-]
%   et_i          [1x1]       Initial ephemeris time for propagation                                [s]
%
% OUTPUT
%   residual_tot  [nx3]       Weighted residuals between predicted and real measurements            [varies]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

n = length(stationNames);
topoFrame = cell(1,n);
for i = 1:n
    topoFrame{i} = [stationNames{i}, '_TOPO']; 
end

% Propagate x to the epochs of the measurements
options = odeset('Reltol',1.e-10,'Abstol',1.e-15);

x_tw = cell(1, size(TW, 2));  

for i = 1:size(TW,2)
    t_span = [et_i TW{i}(1)]; %propagate until time window start
    [~,xx] = ode113(fundynamics,t_span,x0,options); %output is row vector
    xx_tw0 = xx(end,:)';
    [~, x_tw_prop] = ode113(fundynamics, TW{i}(:), xx_tw0, options); % x_tw è x_prop nelle righe dopo
    x_tw{i} = x_tw_prop;  %cell storing of tw propagations
end


meas_pred_tot = cell(1, size(meas_real, 2));  %cell initialization

for s=1:size(TW,2) %auxiliary index for time window selection

    [SMOS_azimuth_pred, SMOS_elevation_pred, SMOS_range_pred, ~] = ...
    antenna_pointing(stationNames{s}, TW{s}, x_tw{s}');

    meas_pred_tot{s} = [SMOS_range_pred',SMOS_azimuth_pred',SMOS_elevation_pred']; %predicted measurament storing in cell
end


% Compute the residual of the measurements and append it to the output
residual_tot = cell(1, size(meas_real, 2));

for s=1:size(meas_pred_tot,2)
    residual = zeros(size(meas_pred_tot{s}));
    for k=1:length(TW{s})
        diff_meas_weighted = W_m * [meas_pred_tot{s}(k,1) - meas_real{s}(k,1), angdiff(meas_pred_tot{s}(k,2), ...
                                    meas_real{s}(k,2)),angdiff(meas_pred_tot{s}(k,3),meas_real{s}(k,3))]';
        residual(k,:) = diff_meas_weighted';
    end
    residual_tot{s} = residual;
end

% residual_tot = [residual_tot{1};residual_tot{2};residual_tot{3}];
residual_tot = vertcat(residual_tot{:});
end

function [range_real, azimuth_real, elevation_real] = measuraments(TW,et0,stationName,satrec,ddpsi,ddeps)
% DESCRIPTION
%   This function computes simulated noisy range, azimuth, and elevation 
%   measurements of a satellite from a ground station during specified 
%   visibility time windows. It performs SGP4 propagation, converts states 
%   from TEME to ECI frames, then to a topocentric frame, and finally derives 
%   measurement values with added Gaussian noise.
%
% INPUT
%   TW          [nx1]     Vector of ephemeris times (ET) during visibility windows               [s]
%   et0         [1x1]     Reference ephemeris time (ET)                                          [s]
%   stationName [char]    Ground station name as registered in SPICE kernel                     [-]
%   satrec      [struct]  SGP4 satellite record structure                                        [-]
%   ddpsi       [1x1]     Nutation correction in longitude                                       [rad]
%   ddeps       [1x1]     Nutation correction in obliquity                                       [rad]
%
% OUTPUT
%   range_real      [nx1]     Simulated noisy range measurements                                 [km]
%   azimuth_real    [nx1]     Simulated noisy azimuth measurements                               [rad]
%   elevation_real  [nx1]     Simulated noisy elevation measurements                             [rad]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

npoints = length(TW);

%Matrix inizialization to store propagated states a measuraments
rteme = zeros(3,npoints);
vteme = zeros(3,npoints);
reci_SGP4 = zeros(3,npoints);
veci_SGP4 = zeros(3,npoints);
sat_range = zeros(npoints,1);
sat_azimuth  = zeros(npoints,1);
sat_elevation = zeros(npoints,1);
% Define station name
topoFrame = [stationName, '_TOPO'];
fprintf('Station Name: %s\n',topoFrame);

for i = 1:npoints

    % SGP4 propagation
    tsince = (TW(i) - et0)/60; %required in minutes
    [~,rteme(:,i),vteme(:,i)] = sgp4(satrec,  tsince); %SGP works on teme and needs conversion
    
    % Compute centuries from 2000-01-01T00:00:00.00 TDT
    ttt = cspice_unitim(TW(i), 'ET', 'TDT')/cspice_jyear()/100;

    % TEME to ECI conversion
    [reci_SGP4(:,i), veci_SGP4(:,i)] = ...
        teme2eci(rteme(:,i), vteme(:,i), [0.0;0.0;0.0],  ttt, ddpsi, ddeps);

    % Rotation matrix from ECI to topocentric frame
    ROT_ECI2TOPO = cspice_sxform('J2000', topoFrame, TW(i));

    % Compute station position in ECI
    rv_station_eci = cspice_spkezr(stationName, TW(i), 'J2000', 'NONE', 'EARTH');

    % Compute station-satellite vector in ECI
    rv_station_sat_eci = [reci_SGP4(:,i); veci_SGP4(:,i)] - rv_station_eci;

    % Convert state into topocentric frame
    rv_station_sat_topo = ROT_ECI2TOPO*rv_station_sat_eci;

    % Compute range, azimuth and elevation using cspice_xfmsta
    rll_station_sat = cspice_xfmsta(rv_station_sat_topo,'RECTANGULAR','LATITUDINAL','EARTH');

    sat_range(i,1)   = rll_station_sat(1);   % [km]
    sat_azimuth(i,1) = rll_station_sat(2);   % [rad]
    sat_elevation(i,1) = rll_station_sat(3); % [rad]
end

sigma_range = 0.01;                %[km]
sigma_azimuth = deg2rad(0.125);    %[rad]
sigma_elevation = deg2rad(0.125);  %[rad]

range_noise_cov = sigma_range^2;
range_real = mvnrnd(sat_range,range_noise_cov);

azimuth_noise_cov = sigma_azimuth^2;
azimuth_real = mvnrnd(sat_azimuth,azimuth_noise_cov);

elevation_noise_cov = sigma_elevation^2;
elevation_real = mvnrnd(sat_elevation,elevation_noise_cov);
end

function [dxdt] = J2_rhs(~, x, mu)
% DESCRIPTION:
% This function describes the right hand side of two body dynamics plus the
% perturbation caused by second zonal spherical harmonic (J2)
%
%
%
% Inputs:
%   t   : [ 1, 1] epoch (unused)
%   x   : [6, 1] cartesian state vector wrt ECI
%   mu  : [ 1, 1] gravitational constant of the body
%
% Outputs:
%   dxdt   : [6,1] RHS, newtonian gravitational acceleration and second zonal spherical harmonic


r=x(1:3);

rnorm   = norm(r);

%Define J2 constant
J2=1082.63*1e-6; %taken from (https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)

%Define Radius of Earth (Average) 
Re = cspice_bodvrd('399','RADII',3);
Re = (Re(1)+Re(2)+Re(3))/3;

%compute acceleration due to J2 in x,y and z direction
J2_acc(1)     = 1.5*((J2*mu*Re^2)/rnorm^4)*(r(1)/rnorm)*(5*(r(3)^2/rnorm^2)-1);%x
J2_acc(2)     = 1.5*((J2*mu*Re^2)/rnorm^4)*(r(2)/rnorm)*(5*(r(3)^2/rnorm^2)-1);%y  
J2_acc(3)    = 1.5*((J2*mu*Re^2)/rnorm^4)*(r(3)/rnorm)*(5*(r(3)^2/rnorm^2)-3);%z

% Compute square distance and distance
dist2 = dot(r, r);
dist = sqrt(dist2);


dxdt=zeros(6,1);%inizialize as column vector
% Position detivative is object's velocity
dxdt(1:3) = x(4:6)';   
% Compute the gravitational acceleration using Newton's law
dxdt(4:6) = - mu * r /(dist*dist2)+J2_acc';

end

function [dxdt] = keplerian_rhs(~, x, GM)
% DESCRIPTION
% This function evaluates the right-hand side of the Newtonian two-body (Keplerian) equations of motion. 
% It returns the derivatives of a 6-dimensional state vector, which includes the object's position and velocity, 
% for propagation in a gravitational field.
%
% INPUT
% t  [1x1]  time epoch (not used in computation)                       [s]
% x  [6x1]  Cartesian state vector: position (1:3) and velocity (4:6)  [km, km/s]
% GM [1x1]  gravitational parameter (G*M) of the primary body          [km^3/s^2]
%
% OUTPUT
% dxdt [6x1] time derivative of the state vector: velocity and acceleration [m/s, m/s^2]
%
% AUTHOR
% Mina Baniamein, A.Y. 2024/25, MATLAB
%

% Initialize right-hand-side
dxdt = zeros(6,1);

% Extract positions
rr = x(1:3);

% Compute square distance and distance
dist2 = dot(rr, rr);
dist = sqrt(dist2);

% Position detivative is object's velocity
dxdt(1:3) = x(4:6);   
% Compute the gravitational acceleration using Newton's law
dxdt(4:6) = - GM * rr /(dist*dist2);

end
