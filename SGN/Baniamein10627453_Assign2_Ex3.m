%Spacecraft Guidance and Navigation
%Assignment #2
%Author: Mina Baniamein

%% Point 1
clear ; close all; clc;
%% Kernel unload
cspice_kclear(); %unload previous kernels section
%% Meta kernel
cspice_furnsh('kernels\assignment2.txt'); 
%% Point 1
clc
%Satellite initial state
r0mci = [4307.844185282820;-1317.980749248651;2109.210101634011];
v0mci = [-0.110997301537882;-0.509392750828585;0.815198807994189];
x0mci = [r0mci;v0mci];
sat_epoch_str_t0 = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [2024,11,18,16,30,0]);
et0 = cspice_str2et(sat_epoch_str_t0);
sat_epoch_str_tf = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [2024,11,18,20,30,0]);
etf = cspice_str2et(sat_epoch_str_tf);
sigma_p = 0.1; %[km]
P0 = diag([10,1,1,0.001,0.001,0.001,0.00001,0.00001]);
mu = cspice_bodvrd('MOON', 'GM', 1);
%Keplerian motion prediction of satellite
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[~, xx_sat] = ode113(@(t,x) keplerian_rhs(t, x, mu), 0:30:etf-et0, x0mci, options);

%Lander position in time
%Define longitude, latitude, and altitude:
lon =  15.0 * cspice_rpd(); % [rad]
lat =   78.0 * cspice_rpd(); % [rad]
alt =   0;  
% Define the time instants
times = et0:30:etf;
% Extract Moons's Radii:
radii_Moon = cspice_bodvrd('MOON','RADII',3);
% Compute Earth's flatness:
f = (radii_Moon(1)-radii_Moon(3))/radii_Moon(1);
% Compute rectangular coordinates from geodetic:
rr_land_rec = cspice_georec(lon, lat, alt, radii_Moon(1), f);
% Initialize storing positions vwctor
rr_lander_eci = zeros(3,length(times));
% Compute rotation matrix to J2000 for every time instant:
for i=1:length(times)
    ROT_EF_J2000 = cspice_pxform('IAU_MOON', 'J2000', times(i));
    % Rotated position vectors:
    rr_lander_eci(:,i) = ROT_EF_J2000*rr_land_rec;
end
% Time window computation
xx_sat = xx_sat';%for easier indexing

%Measurements storing vectors
sat_range = zeros(length(times),1);
sat_azimuth  = zeros(length(times),1);
sat_elevation = zeros(length(times),1);

for i=1:length(times)
    
    % Compute station-lander relative position vector in ECI
    rr_lander_sat_eci = xx_sat(1:3,i) - rr_lander_eci(:,i);
 
    % Convert state into topocentric frame
    [range, azimuth, elevation] = cspice_reclat(rr_lander_sat_eci);

    sat_range(i,1)   = range;       % [km]
    sat_azimuth(i,1) = azimuth;     % [rad]
    sat_elevation(i,1) = elevation; % [rad]

end

i_visibility_lander = sat_elevation > 0;

TW = times(i_visibility_lander); 

if length(TW) == length(times)
    disp('The satellite is visible by the lander during the whole time window');
else
    disp('The satellite is NOT visible by the lander during the whole time window');
end

TW_beginning = cspice_et2utc(TW(1) ,'C', 3);
TW_ending = cspice_et2utc(TW(end) ,'C', 3);
disp(['Visibility window beginning: ', num2str(TW_beginning)]);
disp(['Visibility window ending: ', num2str(TW_ending)]);


%% Point 2
clc
% Topocentric frame of the lander
topoFrame = ['MOONLANDER', '_TOPO']; %exploiting a more precise kernel for measurament simulaton
% Noise
sigma = 0.1; %[km]
range_real = zeros(length(times),1);

for i = 1:length(TW)
    % Rotation matrix from ECI to topocentric frame
    ROT_ECI2TOPO = cspice_sxform('J2000', topoFrame, TW(i));

    % Compute station position in ECI
    rv_lander_eci = cspice_spkezr('MOONLANDER', TW(i), 'J2000', 'NONE', 'MOON');

    % Compute station-satellite vector in ECI
    rv_station_sat_eci = xx_sat(:,i) - rv_lander_eci;

    % Convert state into topocentric frame
    rv_lander_sat_topo = ROT_ECI2TOPO*rv_station_sat_eci;

    % Compute range, azimuth and elevation using cspice_xfmsta
    rll_lander_sat = cspice_xfmsta(rv_lander_sat_topo,'RECTANGULAR','LATITUDINAL','MOON');

    range_real(i,1)   = rll_lander_sat(1) + sigma*randn;   % [km]

end

real_position = xx_sat(1:3,:) + sigma*randn(size(xx_sat(1:3,:)));

% Graphical verification of gaussian distribution of noise applied with
% randn, by analizing high number of samples
num_samples = 10000;
noise_range = sigma * randn(num_samples, 1);
% Histogram of noise
figure;
histogram(noise_range, 100, 'Normalization', 'pdf');
hold on;
x = linspace(min(noise_range), max(noise_range), 100);
plot(x, normpdf(x, 0, sigma), 'r', 'LineWidth', 1);
title('Relative Range Noise');
xlabel('\sigma')
ylabel('f(x|\mu,\sigma)')
%% Point 3
x0 = (mvnrnd(x0mci,P0(1:6,1:6)))'; 
alpha = 0.01;
beta = 2;
R = sigma^2 * eye(3);
error = zeros(1,length(TW) - 1);
sigma_posvel = zeros(2, length(TW) - 1);
P = P0(1:6,1:6);
x = x0; 

for i=1:(length(TW)-1)
    [x_post, P_post] = UKF1(x, TW(i), TW(i+1), P, alpha, beta, mu, R, real_position(:,i+1));
    error(i) = norm(xx_sat(1:3,i+1) - x_post(1:3));
    sigma_posvel(:,i) = [3*sqrt(trace((P_post(1:3,1:3)))); 3*sqrt(trace((P_post(4:6,4:6))))]; 
    x = x_post;
    P = P_post;
end

figure
plot(TW(2:end),error)
title('Error vs time')
xlabel('t')
ylabel('error')

figure
hold on
plot(TW(2:end), sigma_posvel(1,:), 'r');
plot(TW(2:end), sigma_posvel(2,:), 'b');
title('3\sigma vs time')
legend('3\sigma_r','3\sigma_v')
xlabel('t')
ylabel('3\sigma')

%% Point 4
clc
x0 = (mvnrnd([x0mci; lat; lon],P0))'; 
alpha = 0.01;
beta = 2;
R = sigma^2 * eye(4);
error2 = zeros(1,length(TW) - 1);
sigma_posvelcoord = zeros(4, length(TW) - 1);
P = P0;
x = x0;

for i=1:(length(TW)-1)
    [x_post, P_post] = UKF2(x, TW(i), TW(i+1), P, alpha, beta, mu, R, [real_position(:,i+1); range_real(i+1)]);
    error2(i) = norm( [xx_sat(1:3,i+1) - x_post(1:3); lat - x_post(7); lon - x_post(8)]);
    sigma_posvelcoord(1:2,i) = [3*sqrt(trace((P_post(1:3,1:3)))); 3*sqrt(trace((P_post(4:6,4:6))))];
    sigma_posvelcoord(3,i) = 3*sqrt(P_post(7,7));
    sigma_posvelcoord(4,i) = 3*sqrt(P_post(8,8));
    x = x_post;
    P = P_post;
end

figure
plot(TW(2:end),error2)
title('Error vs time')
xlabel('t')
ylabel('error')

figure
hold on;
plot(TW(2:end), sigma_posvelcoord(1,:), 'r');
plot(TW(2:end), sigma_posvelcoord(2,:), 'b');
plot(TW(2:end), sigma_posvelcoord(3,:), 'g');
plot(TW(2:end), sigma_posvelcoord(4,:), 'm--');
title('3\sigma vs time')
legend('3\sigma_r','3\sigma_v','3\sigma_\phi','3\sigma_\lambda')
xlabel('t')
ylabel('3\sigma')
%% Functions
function [x_post, P_post] = UKF2(x0,ti,tf,P0,alpha,beta,mu,R,y)
% DESCRIPTION
% This function implements the Unscented Kalman Filter (UKF) prediction and update steps for a nonlinear system
% whose state vector includes both the satellite state and a fixed ground station location (latitude and longitude). 
% The UKF propagates the satellite states via Keplerian dynamics and updates the state vector and covariance matrix 
% using range observations between the satellite and the fixed ground site.
%
% INPUT
% x0    [8x1]  initial augmented state vector (satellite position [1:3], velocity [4:6], ground station latitude [7], longitude [8]) [-, -, rad, rad]
% ti    [1x1]  initial time value                                       [s]
% tf    [1x1]  final time value                                         [s]
% P0    [8x8]  initial covariance matrix of the augmented state         [-]
% alpha [1x1]  spread of sigma points parameter                         [-]
% beta  [1x1]  prior knowledge of the distribution                      [-]
% mu    [1x1]  gravitational parameter                                  [-]
% R     [4x4]  measurement noise covariance matrix                      [-]
% y     [4x1]  measurement vector (satellite position [1:3], range [4]) [-, km]
%
% OUTPUT
% x_post [8x1]  a posteriori augmented state estimate after measurement update [-]
% P_post [8x8]  a posteriori covariance matrix after measurement update        [-]
%
% AUTHOR
% Mina Baniamein, A.Y. 2024/25, MATLAB

n = length(x0);
lambda = n*alpha^2 - n;
P_scaled = sqrtm((n+lambda)*P0);

%chi points computation
chi = zeros(n,2*n +1); %inizialization of a matrix to store all chi column vectors
chi(:,1) = x0; 
for i=1:n
    chi(:,1+i) = x0 + P_scaled(:,i);
    chi(:,1+n+i) = x0 - P_scaled(:,i);
end

% weights
Wm0 = lambda/(n+lambda);
Wc0 = lambda/(n+lambda) + 1-alpha^2+beta;
Wi  = 1/(2*(n+lambda));

Wm = [Wm0, Wi.*ones(1,2*n)];
Wc = [Wc0, Wi.*ones(1,2*n)];

% Sigma points propagation for a priori mean state and covariance
chi_prop = zeros(size(chi));

options = odeset('reltol', 1e-12, 'abstol', 1e-12);
for i=1:(2*n+1)
    [~, xx] = ode113(@(t,x) keplerian_rhs(t,x,mu), [ti tf], chi(1:6,i), options); %final state is row, need to transpose
    chi_prop(1:6,i) = [xx(end,:)]'; %sigma points
end
chi_prop(7:8,:) = chi(7:8,:); %coordinates stay the same during the propagation, because the station is fixed

x_m_priori = zeros(size(x0));
for i=1:(2*n + 1)
    x_m_priori = x_m_priori + Wm(1,i).*chi_prop(:,i);
end

P_priori = zeros(n);
for i=1:(2*n + 1)
    P_priori  = P_priori + Wc(1,i).*((chi_prop(:,i)-x_m_priori)*( chi_prop(:,i)-x_m_priori )');
end

range = zeros(1,size(chi_prop,2));
alt=0;
times = tf;
for i = 1:size(chi_prop,2)
    lat = chi_prop(7,i);
    lon = chi_prop(8,i);
    range(i) = getrange(lat,lon,alt,times,chi_prop(1:6,i));
end

Y = [chi_prop(1:3,:);
     range];
y_m = zeros(4,1);
for i=1:(2*n + 1)
    y_m = y_m + Wm(1,i).*Y(:,i);
end

Pyy = zeros(4);  %inizialization of covariance
Pxy = zeros(8,4);  %inizialization of crosscovariance
for i=1:(2*n + 1)
    Pyy  = Pyy + Wc(1,i).*((Y(:,i)-y_m)*( Y(:,i)-y_m )');
end
Pyy = Pyy + R;
for i=1:(2*n + 1)
    Pxy  = Pxy + Wc(1,i).*((chi_prop(:,i)-x_m_priori)*( Y(:,i)-y_m )');
end
% Kalman gain
K = Pxy / Pyy;
% Posteriori results
x_post = x_m_priori + K*(y - y_m);
P_post = P_priori - K*(Pyy)*K';
end

function range = getrange(lat,lon,alt,times,xx_sat)
% DESCRIPTION
% This function computes the geometric range (distance) between a ground site and a satellite, 
% for a series of time instants. The ground site is defined by geodetic coordinates (latitude, longitude, altitude) on the Moon, 
% and the satellite positions are provided in inertial coordinates for each time instant.
%
% INPUT
% lat    [1x1] geodetic latitude of the ground site                   [rad]
% lon    [1x1] geodetic longitude of the ground site                  [rad]
% alt    [1x1] altitude above the Moon reference ellipsoid            [km]
% times  [1xN] row vector of time instants (ephemeris time, ET)       [s]
% xx_sat [6xN] state vectors (position and velocity) of satellite in ECI frame for each time instant:
%              position (1:3,:), velocity (4:6,:)                     [km, km/s]
%
% OUTPUT
% range  [Nx1] geometric range between ground site and satellite at each time instant [km]
%
% AUTHOR
% Mina Baniamein, A.Y. 2024/25, MATLAB

% Extract Moons's Radii:
radii_Moon = cspice_bodvrd('MOON','RADII',3);
% Compute Earth's flatness:
f = (radii_Moon(1)-radii_Moon(3))/radii_Moon(1);
% Compute rectangular coordinates from geodetic:
rr_land_rec = cspice_georec(lon, lat, alt, radii_Moon(1), f);
% Initialize storing positions vwctor
rr_lander_eci = zeros(3,length(times));
% Compute rotation matrix to J2000 for every time instant:
for i=1:length(times)
    ROT_EF_J2000 = cspice_pxform('IAU_MOON', 'J2000', times(i));
    % Rotated position vectors:
    rr_lander_eci(:,i) = ROT_EF_J2000*rr_land_rec;
end

%Measurements storing vectors
sat_range = zeros(length(times),1);

for i=1:length(times)
    
    % Compute station-lander relative position vector in ECI
    rr_lander_sat_eci = xx_sat(1:3,i) - rr_lander_eci(:,i);
 
    % Convert state into topocentric frame
    [range, ~, ~] = cspice_reclat(rr_lander_sat_eci);
    sat_range(i) = range;

end
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

function [x_post, P_post] = UKF1(x0,ti,tf,P0,alpha,beta,mu,R,y)
% DESCRIPTION
% This function implements the Unscented Kalman Filter (UKF) prediction and update steps for a nonlinear system, 
% specifically for a system whose dynamics are propagated with the Keplerian right-hand side. It computes the 
% a posteriori state estimate and covariance matrix based on given initial estimates, measurement, and system parameters.
%
% INPUT
% x0   [nx1] vector containing the initial state estimate                 [-]
% ti   [1x1] initial time value                                           [s]
% tf   [1x1] final time value                                             [s]
% P0   [nxn] initial state covariance matrix                              [-]
% alpha[1x1] spread of sigma points parameter                             [-]
% beta [1x1] prior knowledge of the distribution                         [-]
% mu   [1x1] gravitational parameter                                      [-]
% R    [3x3] measurement noise covariance matrix                          [-]
% y    [3x1] measurement vector                                           [-]
%
% OUTPUT
% x_post [nx1] a posteriori state estimate after measurement update        [-]
% P_post [nxn] a posteriori error covariance matrix after measurement update[-]
%
% AUTHOR
% Mina Baniamein, A.Y. 2024/25, MATLAB

n = length(x0);
lambda = n*alpha^2 - n;
P_scaled = sqrtm((n+lambda)*P0);

%chi points computation
chi = zeros(n,2*n +1); %inizialization of a 4x9 matrix to store all chi column vectors
chi(:,1) = x0; 
for i=1:n
    chi(:,1+i) = x0 + P_scaled(:,i);
    chi(:,1+n+i) = x0 - P_scaled(:,i);
end

% weights
Wm0 = lambda/(n+lambda);
Wc0 = lambda/(n+lambda) + 1-alpha^2+beta;
Wi  = 1/(2*(n+lambda));

Wm = [Wm0, Wi.*ones(1,2*n)];
Wc = [Wc0, Wi.*ones(1,2*n)];

% Sigma points propagation for a priori mean state and covariance
chi_prop = zeros(size(chi));

options = odeset('reltol', 1e-12, 'abstol', 1e-12);
for i=1:(2*n+1)
    [~, xx] = ode113(@(t,x) keplerian_rhs(t,x,mu), [ti tf], chi(:,i), options); %final state is row, need to transpose
    chi_prop(:,i) = [xx(end,:)]'; %sigma points
end

x_m_priori = zeros(size(x0));
for i=1:(2*n + 1)
    x_m_priori = x_m_priori + Wm(1,i).*chi_prop(:,i);
end

P_priori = zeros(n);
for i=1:(2*n + 1)
    P_priori  = P_priori + Wc(1,i).*((chi_prop(:,i)-x_m_priori)*( chi_prop(:,i)-x_m_priori )');
end

Y = chi_prop(1:3,:);
y_m = zeros(3,1);
for i=1:(2*n + 1)
    y_m = y_m + Wm(1,i).*Y(:,i);
end

Pyy = zeros(3);  %inizialization of covariance
Pxy = zeros(6,3);  %inizialization of crosscovariance
for i=1:(2*n + 1)
    Pyy  = Pyy + Wc(1,i).*((Y(:,i)-y_m)*( Y(:,i)-y_m )');
end
Pyy = Pyy + R;
for i=1:(2*n + 1)
    Pxy  = Pxy + Wc(1,i).*((chi_prop(:,i)-x_m_priori)*( Y(:,i)-y_m )');
end
% Kalman gain
K = Pxy / Pyy;
% Posteriori results
x_post = x_m_priori + K*(y - y_m);
P_post = P_priori - K*(Pyy)*K';
end
