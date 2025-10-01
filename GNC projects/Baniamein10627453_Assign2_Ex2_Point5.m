% Point 5 of SGN exercise 2.5
clear ; close all; clc;
%% Kernel unload
cspice_kclear(); %unload previous kernels section
%% Meta kernel
cspice_furnsh('kernels\assignment2.txt'); 
%%
clc
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
% Evaluate the TLE
[satrec,rteme,vteme] = sgp4(satrec,  0.0);

ddpsi = -0.115178;
ddeps = -0.007395;

et_i = cspice_str2et(sat_epoch_str);  %t ref= UTC 2024-11-18 18:15:16.064640 
ttt = cspice_unitim(et_i, 'ET', 'TDT')/cspice_jyear()/100; % %number of centuries past J2000 for TEMEtoECI conversion 

% Transform TEME to ECI vectors
ateme = [0;0;0];
[reci_i, veci_i, aeci_i] = teme2eci(rteme, vteme, ateme, ttt, ddpsi, ddeps); %reci and veci are 3x1 column

mu = satrec.mu;

sat_epochf_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [2024,11,28,18,15,16]);
et_f=cspice_str2et(sat_epochf_str); %expressed in seconds past J2000 epoch

% Equally spaced vector, with one minute time-step FOR KOUROU STATION AND
% SVALBARD
npoints = round((et_f-et_i)/60.0)+1;    %minutes past J2000
et_vec = linspace(et_i, et_f, npoints); %minutes sized spacing

%Equally spaced vector, with 30s time-step because sampling time of
%TROLL
npoints_troll = round((et_f-et_i)/30.0)+1; %minutes past J2000
et_vec_troll = linspace(et_i, et_f, npoints_troll); %half minutes sized spacing 

%Propagation of 1 day
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[~, xxf] = ode113(@(t,x) J2_rhs(t, x, mu), et_vec, [reci_i; veci_i], options);
[~, xxf_troll] = ode113(@(t,x) J2_rhs(t, x, mu), et_vec_troll, [reci_i; veci_i], options);

reci_f=xxf(end,1:3)';
reci_f_troll=xxf_troll(end,1:3)'; %check that the final state computed is the same with a refined time grid
veci_f=xxf(end,4:6)';

% Compute antenna angles, satellite range and range-rate
[SMOS_azimuth_kourou, SMOS_elevation_kourou, SMOS_range_kourou, SMOS_range_rate_kourou] = ...
    antenna_pointing('KOUROU', et_vec, [xxf(:,1:3)';xxf(:,4:6)']);

[SMOS_azimuth_troll, SMOS_elevation_troll, SMOS_range_troll, SMOS_range_rate_troll] = ...
    antenna_pointing('TROLL', et_vec_troll, [xxf_troll(:,1:3)';xxf_troll(:,4:6)']);

[SMOS_azimuth_svalbard, SMOS_elevation_svalbard, SMOS_range_svalbard, SMOS_range_rate_svalbard] = ...
    antenna_pointing('SVALBARD', et_vec, [xxf(:,1:3)';xxf(:,4:6)']);

i_visibility_kourou = SMOS_elevation_kourou > deg2rad(6);
i_visibility_troll = SMOS_elevation_troll > deg2rad(0);
i_visibility_svalbard = SMOS_elevation_svalbard > deg2rad(8);

%Korou measuraments computation
TW_kourou = et_vec(i_visibility_kourou);
stationName = 'KOUROU';
[range_real_kourou, azimuth_real_kourou, elevation_real_kourou] = measuraments(TW_kourou,et_i,stationName,satrec,ddpsi,ddeps);
%Troll measuraments computation
TW_troll = et_vec_troll(i_visibility_troll);
stationName = 'TROLL';
[range_real_troll, azimuth_real_troll, elevation_real_troll] = measuraments(TW_troll,et_i,stationName,satrec,ddpsi,ddeps);
%Svalbard measuraments computation
TW_svalbard = et_vec(i_visibility_svalbard);
stationName = 'SVALBARD';
[range_real_svalbard, azimuth_real_svalbard, elevation_real_svalbard] = measuraments(TW_svalbard,et_i,stationName,satrec,ddpsi,ddeps);

sigma_range = 0.01;                %[km]
sigma_azimuth = deg2rad(0.125);    %[rad]
sigma_elevation = deg2rad(0.125);  %[rad]
%Weigths matrix
meas_noise_cov = diag([sigma_range^2;sigma_azimuth^2;sigma_elevation^2]);
W_m = inv(sqrtm(meas_noise_cov));

%%
% KOUROU - TROLL
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
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt');
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
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt');
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
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt');
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

%% Commenti per risultati
%1 giorno: TS
%2 giorni: ballottaggio tra stazioni con Troll
%10 giorni: ballottaggio tra stazioni con Troll
%Molto probabilmente: Troll --> più misure ---> più certezza
%Valutare quindi performance delle singole stazioni per capire quale
%scegliere come prime e quale come backup




















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
    % meas_pred = zeros(length(TW{s}),3);

    % for i=1:length(TW{s}) %for cycle compute predictions and store it in a matrix
    %     rv_station_eciTWi = cspice_spkezr(stationNames{s}, TW{s}(i), 'J2000', 'NONE', 'EARTH'); %station position, 6x1
    %     rv_station_sat_eciTWi = x_tw{s}(i,:)' - rv_station_eciTWi; %relative position at beginning of time window
    %     ROT_ECI2TOPO = cspice_sxform('J2000', topoFrame{s}, TW{s}(i));
    %     rv_station_sat_topoTWi = ROT_ECI2TOPO*rv_station_sat_eciTWi;
    %     xx_kepTWi_LAT = cspice_xfmsta(rv_station_sat_topoTWi,'RECTANGULAR','LATITUDINAL','EARTH'); %theoretical measuraments
    %     xx_kepTWi = xx_kepTWi_LAT(1:3); %theoretical measuraments of range,azimuth and elevation
    %     meas_pred(i,:) = xx_kepTWi';    %measuraments storing
    % end
    % 
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
%


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

