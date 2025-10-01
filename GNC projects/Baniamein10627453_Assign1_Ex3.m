%Spacecraft Guidance and Navigation
%Assignment #1
%Author: Mina Baniamein
%% Exercise 3
clear ; close all; clc; 
%% Point 1
% Spatial density plot
%Data
hi = 800;        %[km]
hf = 1000;       %[km]
di = 0.75;       %[deg]
Re = 6378.1366;  %[km]
mu = 398600.435; %[km^3/s^2]
rho0 = 750 + Re; %[km]
k1 = 10^-5;      %[DU^-1]
k2 = 10^-4;      %[DU^2]
m0 = 1000;       %[kg]
Tmax = 3;        %[N]
Isp = 3120;      %[s]
DU = 7178.1366;  %[km]
MU = m0;         %[kg]
g0 = 9.81;       %[m/s^2]

ri = Re + hi;
rf = Re + hf;
rho_plot = linspace(ri - 100, rf + 100, 400);
q_rho_plot = zeros(1,length(rho_plot));


for i=1:length(rho_plot)
    q_rho_plot(i) = spatialdensity(k1,k2,rho_plot(i),rho0,DU);
end

figure
grid on
plot(rho_plot,q_rho_plot)
xlabel('\rho')
ylabel('q(\rho)')
title('Spatial density plot')

%Initial state
v0 = sqrt(mu/ri);
x0 = [ri,0,0,0,v0,0]';
%Final state
vf = sqrt(mu/rf);
xf = [rf,0,0,0,vf*cosd(di),vf*sind(di)]';  
GM = mu;

%% Point 2
% Adimensionalization
% Run this section one time to adimensionalize all the parameters 
GM_ad = 1;
m0_ad = m0/MU;
TU = sqrt((DU)^3/398600.435);
ri_ad = ri/DU;
rf_ad = rf/DU;
rho0_ad = rho0/DU;
Tmax_ad = Tmax/(MU*(DU*1000)/TU^2);
Isp_ad = Isp/TU;
VU = DU/TU;
v0_ad = v0/VU;
vf_ad = vf/VU;
x0_ad = [ri_ad,0,0,0,v0_ad,0]';
xf_ad = [rf_ad,0,0,0,vf_ad*cosd(di),vf_ad*sind(di)]';  
g0_ad = g0*10^-3*(TU^2)/DU;
%% Point 4
% PMP problem solving
% Check Hamiltonian is constant with the implemented dynamics
lambda_rv0 = -250.*ones(6,1) + 500*rand(6,1);
lambda_m0  = 250*rand;
tf0 = (10*2*pi + 2*rand(1) - 1);
lambdatf_0 = [lambda_rv0; lambda_m0; tf0];
parameters = [Tmax_ad,g0_ad,Isp_ad,k1,k2,rho0_ad]; 
options_ode=odeset('RelTol',1e-12,'AbsTol',1e-13);
[t, xx_ham] = ode113(@(t, xx) dynamics(t,xx,GM_ad,parameters), [0 tf0], [x0_ad;m0_ad;lambdatf_0(1:end-1)],options_ode);

H_ham = zeros(length(t),1);
for i=1:length(t)
    H_ham(i) = hamiltonian(k1,k2,xx_ham(i,:)',rho0_ad,GM_ad,Tmax_ad,Isp_ad,g0_ad,1);
end

plot(t,H_ham)
%% Costate solution search
% Solver initialization
counter=0;
exitflag=0;

while exitflag ~= 1
counter=counter+1;
%Guess generation
lambda_r1 = -250 + 500*rand;
lambda_r2 = -250 + 500*rand;
lambda_r3 = -250 + 500*rand;
lambda_v1 = -250 + 500*rand;
lambda_v2 = -250 + 500*rand;
lambda_v3 = -250 + 500*rand;
lambda_m  = 250*rand;

tf0 = (10*2*pi + 2*rand(1) - 1);
z_guess = [lambda_r1; lambda_r2; lambda_r3; lambda_v1; lambda_v2; lambda_v3; lambda_m; tf0];

parameters_ad = [Tmax_ad,g0_ad,Isp_ad,k1,k2,rho0_ad]; 
options_fsolve=optimoptions('fsolve','Algorithm','levenberg-marquardt','FiniteDifferenceType','central','FunctionTolerance', 1e-8, 'TolX', 1e-8, 'MaxFunctionEvaluations', 1e10, 'MaxIterations', 500);
[sol, ~, exitflag] = fsolve(@(z) zero_function(z, x0_ad,m0_ad, xf_ad,parameters_ad),z_guess,options_fsolve);

    % Extract tf from the solution
    tf_sol = sol(end); 
    % If tf is out of the desired range, reject the solution and continue
    if tf_sol < 18*pi || tf_sol > 22*pi
        exitflag = 0;  % Force another iteration
    end

end

%% Solution testing
% The previous section requires long computation time, in this section a
% solution is saved and tested in order to avoid redoing the computations.
% SOLUTION: -214.9812   %lambdas 
%            -10.3659     | 
%              0.8856     |
%            -10.3929     |
%           -214.6105     |
%           -112.9454     |
%              2.5964     V
%             64.4801    %tf

sol = [-214.9812; -10.3659;  0.8856; -10.3929; -214.6105; -112.9454; 2.5964; 64.4801];
options_ode=odeset('RelTol',1e-13,'AbsTol',1e-13);
[t, xx_sol] = ode113(@(t, xx) dynamics(t,xx,GM_ad,parameters), [0 sol(end)], [x0_ad;m0_ad;sol(1:(end-1))],options_ode);
 
% Final state errors in km and m/s
error = xx_sol(end,1:6)' - xf_ad;
error_r = norm(error(1:3).*DU);
error_v = norm(error(4:6).*VU);
disp(['Final error in position: ',num2str(error_r)]);
disp(['Final error in velocity: ',num2str(error_v)]);
% Fuel final mass in kg
m_f = MU*xx_sol(end,7);
disp(['Final fuel mass: ',num2str(m_f)]);
% Final time in minutes
tf_min = 64.4801*TU/60;
disp(['Final time in minutes: ',num2str(tf_min)]);

%% Plot of 3N solution
figure
hold on
grid on
plot3(xx_sol(:,1).*DU, xx_sol(:,2).*DU, xx_sol(:,3).*DU, 'm')
plot3(xx_sol(1,1).*DU, xx_sol(1,2).*DU, xx_sol(1,3).*DU, 'ko', 'MarkerSize', 2,'MarkerFaceColor', 'r')
plot3(xx_sol(end,1).*DU, xx_sol(end,2).*DU, xx_sol(end,3).*DU, 'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'g')
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')
legend('Trajectory','Departure point','Arrival Point'),
title('3N fuel optimal trajectory')

%% Hamiltonian graphical chek of the solution
H_sol = zeros(size(xx_sol,1),1);
options_ode=odeset('RelTol',1e-12,'AbsTol',1e-13);

for i=1:size(xx_sol,1)
    H_sol(i) = hamiltonian(k1,k2,xx_sol(i,:)',rho0_ad,GM_ad,Tmax_ad,Isp_ad,g0_ad,1);
end

figure
plot(t,H_sol)
ylim([-1,1])
title('Hamiltonian value evolution in time')
xlabel('t [-]')
ylabel('H')

%% NTW plot
% Transform alpha into NTW frame
alpha = zeros(3, length(xx_sol)); % Allocate alpha array
state_NTW = zeros(6, length(xx_sol)); % Allocate transformed state array
rot = zeros(3,3,length(xx_sol));
for i = 1:length(xx_sol)
    rot(:,:,i) = ECI2NTW(xx_sol(i,1:3)', xx_sol(i,4:6)');
    alpha(:,i) = rot(:,:,i) * (-xx_sol(i,11:13)' / norm(xx_sol(i,11:13)));
end

% Plot alpha components over time
figure; plot(t, alpha(1,:)); xlim([0 sol(end)]); title('$\alpha_N$ evolution in time','Interpreter', 'latex'); xlabel('Time [s]'); ylabel('\alpha_N [-]');
figure; plot(t, alpha(2,:)); xlim([0 sol(end)]); title('$\alpha_T$ evolution in time','Interpreter', 'latex'); xlabel('Time [s]'); ylabel('\alpha_T [-]');
figure; plot(t, alpha(3,:)); xlim([0 sol(end)]); title('$\alpha_W$ evolution in time','Interpreter', 'latex'); xlabel('Time [s]'); ylabel('\alpha_W [-]');


%% Point 5 
% Numerical continuation
clc
Tmax_f_ad = 2.860 * TU^2 / (MU * DU * 1e3);
T_span = linspace(Tmax_ad, Tmax_f_ad, 100);

% Solver initialization
options_fsolve = optimoptions('fsolve', 'Display', 'iter', ...
    'FunctionTolerance', 1e-6, 'OptimalityTolerance', 1e-6, 'MaxIterations', 100); 

%Solution of Tmax = 3 N serves as first guess for numerical continuation
guess = [-214.9812; -10.3659; 0.8856; -10.3929; -214.6105; -112.9454; 2.5964; 64.4801];

for i = 1:length(T_span)
    parameters(1) = T_span(i); % Update thrust in parameters
    [sol, ~, exitflag] = fsolve(@(z) zero_function(z, x0_ad,m0_ad, xf_ad,parameters),guess,options_fsolve);
    guess = sol; % Use the solution as the next initial guess
end

%% Precompputed solution
% SOLUTION:  -425.9501   %lambdas
%             -11.2531   |
%               1.2846   |
%             -11.4431   |
%            -425.5786   |
%            -547.3406   |
%              10.4493   V
%              64.2507   %tf
clc 
Tmax_f_ad = 2.860*TU^2/(MU*DU* 1e3);
sol = [-425.9501, -11.2531, 1.2846, -11.4431, -425.5786 ,-547.3406, 10.4493, 64.2507]';
parameters(1)  = Tmax_f_ad;
% Extract updated results from final solution
xx0 = [x0_ad; m0_ad; sol(1:end-1)];
lambdas_cont = sol(1:end-1);          % Updated initial costate for mass
t_f_ad_cont = sol(end);               % Updated time of flight [adimensional]
t_f_cont = t_f_ad_cont*TU/60;         % Time of flight [minutes]

% Solve ODE with updated conditions
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[t_cont, x_cont] = ode113(@(t, xx) dynamics(t,xx,GM_ad,parameters), [0, t_f_ad_cont], xx0, options);

% Final state errors in km and m/s
error_cont = x_cont(end,1:6)' - xf_ad;
error_r_cont = norm(error_cont(1:3).*DU);
disp(['Error in position: ', num2str(error_r_cont)]);
error_v_cont = norm(error_cont(4:6).*VU);
disp(['Error in velocity: ', num2str(error_v_cont)]);
% Fuel final mass in kg
m_f_cont = MU*x_cont(end,7);
disp(['Final mass: ',num2str(m_f_cont)]);
% Final time in minutes
tf_min_cont = 64.2507*TU/60;
disp(['Final time in minutes: ',num2str(tf_min_cont)]);

%% Plot of continuated solution
figure
hold on
grid on
plot3(x_cont(:,1).*DU, x_cont(:,2).*DU, x_cont(:,3).*DU, 'm')
plot3(x_cont(1,1).*DU, x_cont(1,2).*DU, x_cont(1,3).*DU, 'ko', 'MarkerSize', 3,'MarkerFaceColor', 'r')
plot3(x_cont(end,1).*DU, x_cont(end,2).*DU, x_cont(end,3).*DU, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'b')
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')
legend('Trajectory','Departure point','Arrival Point')
title('2.860 N fuel optimal trajectory')

%% Plot Hamiltonian of continuated
H_ham_cont = zeros(length(t_cont),1);
for i=1:length(t_cont)
    H_ham_cont(i) = hamiltonian(k1,k2,x_cont(i,:)',rho0_ad,GM_ad,Tmax_ad,Isp_ad,g0_ad,1);
end

figure
plot(t_cont,H_ham_cont)
ylim([-1,1])
xlabel('t [-]')
ylabel('H')
title('2.860 N solution H evolution in time')


%%
%Local functions

function F=zero_function(z, x0, m0, xf_target,parameters)
% DESCRIPTION
%   This function computes the zero vector of boundary conditions and transversality
%   constraints for a low-thrust optimal control problem, based on the final state 
%   error, the final costate condition, and the Hamiltonian transversality condition.
%
% INPUT
%   z            [8x1]     vector of unknowns, including initial costates [7x1] and adimensionalized final time [1x1]                 [-]
%   x0           [6x1]     initial adimensionalized state vector, containing position and velocity components                         [-]
%   m0           [1x1]     initial adimensionalized spacecraft mass                                                                   [-]
%   xf_target    [6x1]     target adimensionalized final state vector, containing position and velocity components                    [-]
%   parameters   [6x1]     vector containing problem parameters: [Tmax ([-]), g0 ([-]), Isp ([-]), k1 ([-]), k2 ([-]), rho0 ([-])]    [-]
%
% OUTPUT
%   F            [8x1]     zero vector of boundary conditions and transversality constraints                                          [-]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

GM = 1;
Tmax = parameters(1);
g0   = parameters(2);
Isp  = parameters(3);
k1   = parameters(4);
k2   = parameters(5);
rho0 = parameters(6);

%z is [lambda0;tf_ad] %unknown of problem
lambda0 = z(1:7);
tf=z(8);

%final state target
rr_f_target = xf_target(1:3);
vv_f_target = xf_target(4:6);


% Propagation of dynamic of state and co-state
%initial condition of state and co-state [rv_earth_ad; m0_ad, lambda_0]
y0=[x0;m0;lambda0];

options_ode=odeset('RelTol',1e-12,'AbsTol',1e-13);
[~,y]=ode113(@(t,y)dynamics(t,y,GM,parameters),[0 tf],y0,options_ode);

%extract final state and co-state
yf=y(end,1:14); %row vector
yf=yf';         %bring back to column vector

%extract useful parameters from state and co-state
rr_f=yf(1:3);
vv_f=yf(4:6);
m_f=yf(7);
r_f=norm(rr_f);


lambda_rr_f=yf(8:10);
lambda_vv_f=yf(11:13);
lambda_m_f=yf(14);
% lambda_f=[lambda_rr_f;lambda_vv_f;lambda_m_f];
% 
% lambda_rr_vv_f=yf(8:13);

% Compute Hamiltonian at final time
%compute derivative of state [r;v;m]
x_dot = [vv_f;
        -GM/(r_f^3)*rr_f-Tmax/m_f*lambda_vv_f/norm(lambda_vv_f);
        -Tmax/(Isp*g0)];

q = spatialdensity(k1,k2,r_f,rho0,1);
H_f=q + dot([lambda_rr_f;lambda_vv_f;lambda_m_f],x_dot);

%Zero function definition
F(1:3)= rr_f-rr_f_target;  %match final position (adimentional)
F(4:6)= vv_f-vv_f_target;  %match final velocity (adimentional)
F(7)  = lambda_m_f;       %lambda_m(tf)=0
F(8)  = H_f;                 %Transversality condition

end

function q = spatialdensity(k1,k2,rho,rho0,DU)
% DESCRIPTION
%   This function calculates the spatial density at a given distance based on
%   a model parameterized by constants k1, k2, a reference distance rho0, and a
%   normalization factor DU. Set DU to 1 for adimensional case
%
% INPUT
%   k1   [1x1]     first model constant coefficient                                                                                [-]
%   k2   [1x1]     second model constant coefficient                                                                                     [-]
%   rho  [1x1]     distance at which to evaluate the spatial density                                                           [L]
%   rho0 [1x1]     reference distance                                                                                         [L]
%   DU   [1x1]     distance unit for normalization (set to 1 if adimensional)                                                   [L or -]
%
% OUTPUT
%   q    [1x1]     computed spatial density value                                                                              [-]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

q = k1 / (k2 + ( (rho-rho0)/DU )^2);

end

function H = hamiltonian(k1,k2,xx,rho0,GM,Tmax,Isp,g0,DU)
% DESCRIPTION
%   This function computes the Hamiltonian value for the low-thrust optimal control
%   problem. Set DU to 1 for adimensional case
%
% INPUT
%   k1     [1x1]     spatial density model constant coefficient                                                              [-]
%   k2     [1x1]     spatial density model constant offset                                                                   [-]
%   xx     [14x1]    combined state and costate vector: 
%                    position (1:3) [L or -], velocity (4:6) [L/T or -], mass (7) [-],
%                    costates lambda_r (8:10) [-], lambda_v (11:13) [-], lambda_m (14) [-]
%   rho0   [1x1]     reference distance for spatial density                                                                   [L or -]
%   GM     [1x1]     gravitational parameter                                                                                  [L^3/T^2 or -]
%   Tmax   [1x1]     maximum thrust                                                                                           [-]
%   Isp    [1x1]     specific impulse                                                                                        [-]
%   g0     [1x1]     standard gravity                                                                                         [-]
%   DU     [1x1]     distance unit for normalization (set to 1 if adimensional)                                                [L or -]
%
% OUTPUT
%   H      [1x1]     scalar value of the Hamiltonian                                                                          [-]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

q = spatialdensity(k1,k2,norm(xx(1:3)),rho0,DU);

H = q + dot(xx(8:10),xx(4:6)) - (GM/(norm(xx(1:3))^3))*dot(xx(1:3),xx(11:13)) + ...
    (Tmax/(Isp*g0))*(-norm(xx(11:13))*Isp*g0/xx(7) - xx(14));
end

function dydt=dynamics(~,y,GM,parameters)
% DESCRIPTION
%   This function computes the time derivatives of the state and costate vectors
%   for the low-thrust optimal control problem dynamics, including spacecraft motion,
%   mass consumption, and costate evolution.
%
% INPUT
%   ~           [1x1]     time variable (not used explicitly)                                                   [-]
%   y           [14x1]    combined state and costate vector:
%                         position (1:3) [L or -], velocity (4:6) [L/T or -], mass (7) [-],
%                         costates lambda_r (8:10) [-], lambda_v (11:13) [-], lambda_m (14) [-]
%   GM          [1x1]     gravitational parameter                                                               [km^3/s^2]
%   parameters  [6x1]     vector containing problem parameters:
%                         [Tmax ([-]), g0 ([-]), Isp ([-]), k1 ([-]), k2 ([-]), rho0 ([-])]
%
% OUTPUT
%   dydt       [14x1]     time derivatives of state and costate vector                                          [-]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

%Constant values
Tmax = parameters(1);
g0   = parameters(2);
Isp  = parameters(3);
k1   = parameters(4);
k2   = parameters(5);
rho0 = parameters(6);

%extract components
rr=y(1:3);
vv=y(4:6);
m=y(7);

lambda_rr=y(8:10);
lambda_vv=y(11:13);
%lambda_m=y(14);

%compute useful quantities
r=norm(rr);
%v=norm(vv);

dydt=zeros(14,1);   %initialize vector

%dx_dt
dydt(1:3)=vv;
dydt(4:6)=-GM/r^3*rr - Tmax/m *lambda_vv/norm(lambda_vv);
dydt(7)=-Tmax/(Isp*g0);

%dqdr AD
dqdr = 2*k1*(r - rho0)/(( (r - rho0)^2 + k2)^2 )*rr/r;

%dlambda_dt
dydt(8:10)=-3*GM/r^5*(rr'*lambda_vv)*rr+GM/r^3*lambda_vv + dqdr;  
dydt(11:13)=-lambda_rr;
dydt(14)=-norm(lambda_vv)*Tmax/(m^2);

end

function rot = ECI2NTW(r, v) 
% DESCRIPTION
% Computes the rotation matrix from ECI to NTW (Normal-Tangential-binormal) frame
%
% INPUT:
%   r   [3x1]  Position vector in ECI frame
%   v   [3x1]  Velocity vector in ECI frame
%
% OUTPUT:
%   rot [3x3]  Rotation matrix from ECI to NTW frame
%
% AUTHOR:
%  Mina Baniamein, 2024/25, MATLAB

    % Compute angular momentum vector
    h = cross(r, v);

    % Compute unit vectors in NTW frame
    N =   h / norm(h);              % Normal direction
    T =   v / norm(v);              % Tangential direction
    W = - r / norm(r);              % Binormal direction
    % Construct rotation matrix
    rot = [N, T, W]';
end
