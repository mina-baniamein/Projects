%Spacecraft Guidance and Navigation
%Assignment #1
%Author: Mina Baniamein

%% Exercise 2
clear ; close all; clc;
%% 
cspice_kclear(); %unload kernel section
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));
%%
cspice_furnsh('kernels\naif0012.tls'); 
cspice_furnsh('kernels\de425s.bsp');   
cspice_furnsh('kernels\gm_de432.tpc');
cspice_furnsh('kernels\pck00010.tpc');
fprintf('Number of LSK  kernels: %d\n', cspice_ktotal('lsk'))
fprintf('Number of SPK  kernels: %d\n', cspice_ktotal('spk'))
fprintf('Number of PCK  kernels: %d\n', cspice_ktotal('pck'))
fprintf('Number of CK   kernels: %d\n', cspice_ktotal('ck'))
fprintf('Number of TEXT kernels: %d\n', cspice_ktotal('TEXT'))
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'))
%% Point 1
% Data
alpha = 0.2*pi;
beta = 1.41;
delta = 4;
t_i = 2;
R_e_vector = cspice_bodvrd('EARTH', 'RADII', 3);
Re = R_e_vector(1); %planar case
R_m_vector = cspice_bodvrd('MOON', 'RADII', 3);
R_moon = R_m_vector(1); %planar case
hi = 167;
hf = 100;
DU = 3.84405000*10^5;
TU = 4.34256461;
VU = 1.02454018 ; 
gm_moon = cspice_bodvrd('MOON', 'GM', 1);
gm_e = cspice_bodvrd('EARTH', 'GM', 1);
mu = gm_moon/(gm_moon + gm_e);
ri = (Re + hi)/DU;
rf = (R_moon + hf)/DU;
r0 = ri;
v0 = beta*sqrt((1-mu)/r0);
ws= -9.25195985*1e-1;

% Guess computation
x0 = r0*cos(alpha) - mu;
y0 = r0*sin(alpha);
vx0 = -(v0 -r0)*sin(alpha);
vy0 = (v0 - r0)*cos(alpha);

xx0 = [x0, y0, vx0, vy0]';

% Guess propagation in the rotating frame
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt, xx] = ode113(@(t,x) PBRFBP(t,x,mu), [2 6], xx0, options);

figure('Name','Rotating Frame') 
plot(xx(:,1),xx(:,2))
hold on 
grid on
axis equal
earth_x = -mu; % x coordinate of Earth 
moon_x = 1-mu; % x coordinate of Moon
plot(earth_x, 0, 'ko', 'MarkerSize', 5,'MarkerFaceColor', '#0019D1'); % Earth
plot(moon_x, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', '#7C8886'); % Moon
circle(earth_x,0,ri); % Initial orbit
circle(moon_x,0,rf); % Final orbit
xlabel('x [-]');
ylabel('y [-]');
title('Rotating frame propagation')
legend('Propagated trajectory','Earth','Moon', 'Initial orbit', 'Final orbit')

% Guess propagation in ECI frame
n = size(xx,1);
xx_ECI = zeros(n,4);

for i=1:n
    ECI = rot_to_ECI(xx(i,:), mu, tt(i));
    xx_ECI(i,1) = ECI(1);
    xx_ECI(i,2) = ECI(2);
    xx_ECI(i,3) = ECI(3);
    xx_ECI(i,4) = ECI(4);
end

earth_ECI = rot_to_ECI([earth_x,0,0,0], mu, tt(1));
moon_ECI  = rot_to_ECI([moon_x,0,0,0], mu, tt(end));
figure('Name','ECI') 
hold on 
grid on
axis equal
plot(xx_ECI(:,1),xx_ECI(:,2));
plot(earth_ECI(1), 0, 'ko', 'MarkerSize', 5,'MarkerFaceColor', '#0019D1'); % Earth
plot(moon_ECI(1), 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', '#7C8886'); % Moon
circle(earth_ECI(1),0,ri); % Initial orbit
circle(moon_ECI(1),0,rf); % Final orbit
xlabel('x [-]');
ylabel('y [-]');
title('ECI frame propagation');
legend('Propagated trajectory','Earth','Moon', 'Initial orbit', 'Final orbit');

% Check if frame change is correct
err = zeros(n,1);
for i=1:n
    err(i) = abs(norm([xx(i,1)+mu,xx(i,2)]) - norm([xx_ECI(i,1), xx_ECI(i,2)]));
end

% Error plot
figure('Name','err')
grid on
plot(linspace(2,6,n), err)
xlabel('t [-]');
ylabel('error [-]');
title('Norm error between frames')

%% Point 2.a
% NLP without gradients

% Initial condition
xx0 = [xx0; 2; 6];

% Fmincon solution research
options=optimset('Algorithm','active-set','LargeScale','on','Display','iter','TolCon',1e-10);  
[xx0_optimal, delta_v] = fmincon(@deltavfun,xx0,[],[],[],[],[],[],@nonlcon,options); %NLP solving

% Resulting initial condition delta v
disp('');
disp(['No gradients solution: r = (' , num2str(xx0_optimal(1)),',', num2str(xx0_optimal(2)), ') v = (', num2str(xx0_optimal(3)),',',num2str(xx0_optimal(4)),')'])
disp(['No gradients dv: ', num2str(delta_v)]);

% Plot of the resulting propagated initial condition
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[~, xx_optimal] = ode113(@(t,x) PBRFBP(t,x,mu), [xx0_optimal(5) xx0_optimal(6)], xx0_optimal(1:4), options); 
figure() 
plot(xx_optimal(:,1),xx_optimal(:,2))
hold on
axis equal
plot(earth_x, 0, 'ko', 'MarkerSize', 5,'MarkerFaceColor', '#0019D1'); % Earth
plot(moon_x, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', '#7C8886'); % Moon
circle(earth_x,0,ri); % Initial orbit
circle(moon_x,0,rf); % Final orbit
xlabel('x [-]');
ylabel('y [-]');
title('Simple shooting trajectory (no gradients)')
legend('Propagated trajectory','Earth','Moon', 'Initial orbit', 'Final orbit')

%% Point 2.b
% NLP with gradients
clc 
% Solver initialization
options=optimset('Algorithm','active-set','LargeScale','on','Display','iter','TolCon',1e-10, ...
    'MaxIter',100,'MaxFunEvals',1000,'GradObj','on','GradConstr','on');

% Fmincon
[xx0_optimal_grad, dv_grad]=fmincon(@deltavfun_grad, xx0, [],[],[],[],[],[],@nonlcon_grad, options);

% Resulting initial condition and delta v
disp('');
disp(['No gradients solution: r = (' , num2str(xx0_optimal_grad(1)),',', num2str(xx0_optimal_grad(2)), ...
    ') v = (', num2str(xx0_optimal_grad(3)),',',num2str(xx0_optimal_grad(4)),')',' ti = ',num2str(xx0_optimal_grad(5)),' tf = ', num2str(xx0_optimal_grad(6))])
disp(['No gradients dv: ', num2str(delta_v)]);

% Plot of the resulting propagated initial condition
figure() 
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[~, xx_optimal_grad] = ode113(@(t,x) PBRFBP(t,x,mu), [xx0_optimal_grad(5) xx0_optimal_grad(6)], xx0_optimal_grad(1:4), options); 
plot(xx_optimal_grad(:,1),xx_optimal_grad(:,2))
hold on
axis equal
plot(earth_x, 0, 'ko', 'MarkerSize', 5,'MarkerFaceColor', '#0019D1'); % Earth
plot(moon_x, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', '#7C8886'); % Moon
circle(earth_x,0,ri); % Initial orbit
circle(moon_x,0,rf); % Final orbit
xlabel('x [-]');
ylabel('y [-]');
title('Simple shooting trajectory (with gradients)')
legend('Propagated trajectory','Earth','Moon', 'Initial orbit', 'Final orbit')
%% Point 3
%  Mutiple shooting NLP with gradients
clc 
% Initial condition
times = linspace(2,6,4);
N = 4; % number of shots
Xopt1=[xx0(1:2); xx0(3:4)];
Xopt0 = zeros(4*N+2,1);
Xopt0(1:4) = Xopt1;
index = 1;
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
for i = 1 : N-1 % Nodes cycle generation
    [~, xx_sol] = ode113(@(t,xx) PBRFBP(t,xx,mu), [times(i), times(i+1)], Xopt0(4*i-3 : 4*i), options);
    Xopt0(4*i+1:4*i+4) = xx_sol(end,:)';
end
Xopt0(end-1) = 2;
Xopt0(end) = 6; 

% Solver intializaion
options=optimset('Algorithm','active-set','LargeScale','on','Display','iter','TolCon',1e-7, ...
    'MaxIter',500,'MaxFunEvals',5000,'GradObj','on','GradConstr','on'); % necessaria una'TolCon'di 1e-7 per farlo convergere in poco tempo

% Fmincon
[xx0_optimal_MS,dv_MS] = fmincon(@deltavfun_MS_grad,Xopt0,[],[],[],[],[],[],@nonlcon_MS_grad,options); 

% Resulting initial condition and delta v
disp('');
disp(['No gradients solution: r = (' , num2str(xx0_optimal_MS(1)),',', num2str(xx0_optimal_MS(2)), ...
    ') v = (', num2str(xx0_optimal_MS(3)),',',num2str(xx0_optimal_MS(4)),')',' ti = ',num2str(xx0_optimal_MS(end-1)),' tf = ', num2str(xx0_optimal_MS(end))])
disp(['No gradients dv: ', num2str(dv_MS)]);

% Plot of the resulting propagated initial condition
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[~, xx_optimal_MS] = ode113(@(t,x) PBRFBP(t,x,mu), [xx0_optimal_MS(end-1) xx0_optimal_MS(end)], xx0_optimal_MS(1:4), options); 

figure() 
plot(xx_optimal_MS(:,1),xx_optimal_MS(:,2))
hold on
axis equal
plot(xx0_optimal_MS(1),xx0_optimal_MS(2),'r.')
plot(xx0_optimal_MS(5),xx0_optimal_MS(6),'r.')
plot(xx0_optimal_MS(9),xx0_optimal_MS(10),'r.')
plot(xx0_optimal_MS(13),xx0_optimal_MS(14),'r.')
plot(earth_x, 0, 'ko', 'MarkerSize', 5,'MarkerFaceColor', '#0019D1'); % Earth
plot(moon_x, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', '#7C8886'); % Moon
circle(earth_x,0,ri,'m'); % Initial orbit
circle(moon_x,0,rf,'g'); % Final orbit
xlabel('x [-]');
ylabel('y [-]');
title('Multiple Shooting trajectory')
legend('Propagated trajectory','Nodes','','', '','Earth','Moon','Initial orbit','Final orbit')
%% Point 4
clc
% Zero-finding problem for departure date

% Data conversion
ti_sol = xx0_optimal_grad(5)*TU*24*3600; 
tf_sol = xx0_optimal_grad(6)*TU*24*3600; 
dt = tf_sol-ti_sol;
theta_i =abs( ws * xx0_optimal_grad(5));

%% Zero function plot
lower_epoch = cspice_str2et('2024-Sep-28 00:00:00.0000 TDB');
upper_epoch = cspice_str2et('2024-Oct-17 00:00:00.0000 TDB');
a = lower_epoch;
b = upper_epoch;
es = linspace(a,b,10000);
fs_storing = zeros(size(es));


for i=1:length(es)
    fs = epoch_search(es(i)) - theta_i;
    fs_storing(i) = fs;
end 

figure
plot(es, fs_storing); xlabel('Epoch'); ylabel('epoch\_search');

%% Zero-finding problem solving
epoch_0 = '2024-Sep-28 00:00:00.0000 TDB';
epoch_0 = cspice_str2et(epoch_0);
lower_epoch = cspice_str2et('2024-Sep-28 00:00:00.0000 TDB');
upper_epoch = cspice_str2et('2024-Oct-17 00:00:00.0000 TDB');
option_fzero = optimset('Display','iter','TolFun',1e-10,'TolX',1e-10);


eti = fzero(@(t) (epoch_search(t) - theta_i), [lower_epoch upper_epoch]);

etf = eti + dt;
epoch_i_str = cspice_et2utc(eti,'C', 3);
epoch_f_str = cspice_et2utc(etf,'C', 3);
disp(['Departure date: ',epoch_i_str])
disp(' ')
disp(['Arrival date: ',epoch_f_str])

%% N-body propagation
% Data initialization
xx_sol_Nbody = zeros(1,6);
% t_conv = 7.820277050997294e+08; % '2024 OCT 12 17:53:55.917'
% str = '2024 OCT 12 17:50:55.917'; 
% eti = cspice_str2et(str);
% 
% % epoch_i_str = cspice_et2utc(t_conv,'C', 3);

% t_conv = eti;
% etf = eti + dt;

t_conv = xx0_optimal_grad(5);
% etf = t_conv + dt;

%%
xx_sol_ECI = rot_to_ECI(xx0_optimal_grad(1:4), mu,t_conv);
xx_sol_Nbody(1:3) = [xx_sol_ECI(1:2),0]*DU;
xx_sol_Nbody(4:6) = [xx_sol_ECI(3:4),0]*VU;
xx_sol_Nbody = xx_sol_Nbody';

% N-body propagation
bodies = {'Sun';'Mercury';'Venus';'Earth';'Moon';'Mars Barycenter';'Jupiter Barycenter';
            'Saturn Barycenter';'Uranus Barycenter';'Neptune Barycenter';'Pluto Barycenter'};

frame = 'ECLIPJ2000'; center = 'Earth';

options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt, xx] = ode78(@(tt,xx) nbody(tt,xx,bodies,frame,center), [eti etf], xx_sol_Nbody, options);



% Plot of the resulting trajectories
figure
axis equal
title ({'Propagation of a transfer orbit Earth-Moon with n-body propagator'},{['@',center,' : ',frame]});
hold on
plot3(xx(:,1), xx(:,2), xx(:,3),'Color','b');
rr_Moon = zeros(3,length(tt));
rr_Earth = zeros(3,length(tt));
for i = 1:length(tt)
    rv_Moon_t = cspice_spkezr('Moon',tt(i),frame,'NONE','Earth');
    rr_Moon(:,i) = rv_Moon_t(1:3);
    rv_Earth_t = cspice_spkezr('Earth',tt(i),frame,'NONE',center);
    rr_Earth(:,i) = rv_Moon_t(1:3);
end
theta = linspace(0,2*pi,1000); 
plot3(0,0,0,'o','Marker', 'o', 'MarkerFaceColor','#33C38A','MarkerEdgeColor','k','MarkerSize',4);
plot3(rr_Moon(1,end), rr_Moon(2,end), rr_Moon(3,end),'o','Marker', 'o', 'MarkerFaceColor','#B7E1FA','MarkerEdgeColor','k','MarkerSize',4);
plot3(rr_Moon(1,:), rr_Moon(2,:), rr_Moon(3,:),'Color','r');
xlabel('x [km]',FontSize=14); ylabel('y [km]',FontSize=14); zlabel('z [km]',FontSize=14);
legend('Transfer Orbit','Earth','Moon','Moon Orbit around Earth','Location','best') ;


%% Functions

function thetaECI = epoch_search(t)
% DESCRIPTION
%   This function computes the angle between the Sun and the Moon, projected 
%   onto the ECLIPJ2000 equatorial plane, as seen from the Earth-Moon Barycenter (EMB) 
%   at a given ephemeris time.
%
% INPUT
%   t           [1x1]       ephemeris time (ET SPICE), seconds past J2000 (TDB)         [s]
%
% OUTPUT
%   thetaECI    [1x1]       angle between the Sun and the Moon in the ECLIPJ2000 
%                           equatorial plane, positive counterclockwise from Sun to Moon [rad]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

    rv_s = cspice_spkezr('MOON', t, 'ECLIPJ2000', 'NONE', 'EMB');
    rv_s(3) = 0;
    r_s = rv_s(1:3)/norm(rv_s(1:3));

    rv_m = cspice_spkezr('SUN', t, 'ECLIPJ2000', 'NONE', 'EMB');
    rv_m(3) = 0;
    r_m = rv_m(1:3)/norm(rv_m(1:3));

    thetaECI = acos(dot(r_m,r_s));

    crossVect = cross(r_m,r_s);

    if crossVect(3) < 0
        thetaECI = - thetaECI;
    end

end

function [dxdt] = PBRFBP(t,xx, mu)
% DESCRIPTION
%   This function contains the RHS for solving the PBRFBP differential
%   equations and computing the state transition matrix
% INPUT
%   ~   [1x1]       placeholder for time input (unused)                                  [s]
%   xx  [4x1]     state vector containing adimensionalized position and velocity vectors [-]
%   mu  [1x1]     adimensionalized mass parameter                                        [-]
% OUTPUT
%   dx  [6x1]     vector containing derivative of the state vector o                     [-] 
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB
    
    % Extract variables
    x  = xx(1);
    y  = xx(2);
    vx = xx(3);
    vy = xx(4);

    % Compute distances from bodies 1 and 2
    r1 = sqrt((x + mu)^2 + y^2);
    r2 = sqrt((x + mu - 1)^2 + y^2);
    
    rho=3.88811143*1e2;
    ms=3.28900541*1e5;
    ws=-9.25195985*1e-1;
    
    dxdt=zeros(4,1);
    dxdt(1) = vx;
    dxdt(2) = vy;
    dxdt(3) = 2*vy+x-(1-mu)*(x+mu)/r1^3-mu*(x+mu-1)/r2^3;
    dxdt(4) = -2*vx+y-(1-mu)*y/r1^3-mu*y/r2^3;

    %PCRTBP perturbation
    r3=sqrt((x-rho*cos(ws*t))^2 + (y-rho*sin(ws*t))^2);
    dxdt(3)=dxdt(3)-ms*(x-rho*cos(ws*t))/r3^3-ms*cos(ws*t)/rho^2;
    dxdt(4)=dxdt(4)-ms*(y-rho*sin(ws*t))/r3^3-ms*sin(ws*t)/rho^2;
end

function [dxdt] = PBRFBP_STM(t,xx, mu)
% DESCRIPTION
%   This function contains the RHS for solving the PBRFBP differential
%   equations and computing the state transition matrix
% INPUT
%   xx  [4x1]     state vector containing adimensionalized position and velocity vectors [-]
%   mu  [1x1]     adimensionalized mass parameter                                        [-]
% OUTPUT
%   dx  [20x1]     vector containing derivative of the state vector o                    [-] 
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB
    % Extract variables
    x  = xx(1);
    y  = xx(2);
    vx = xx(3);
    vy = xx(4);

    PHI = reshape(xx(5:end), 4, 4);

    % Compute distances from bodies 1 and 2
    r1 = sqrt((x + mu)^2 + y^2);
    r2 = sqrt((x + mu - 1)^2 + y^2);
    
    rho=3.88811143*1e2;
    ms=3.28900541*1e5;
    ws=-9.25195985*1e-1;
    
    dxdt=zeros(20,1);
    dxdt(1) = vx;
    dxdt(2) = vy;
    dxdt(3) = 2*vy+x-(1-mu)*(x+mu)/r1^3-mu*(x+mu-1)/r2^3;
    dxdt(4) = -2*vx+y-(1-mu)*y/r1^3-mu*y/r2^3;

    %PCRTBP perturbation
    r3=sqrt((x-rho*cos(ws*t))^2 + (y-rho*sin(ws*t))^2);
    dxdt(3)=dxdt(3)-ms*(x-rho*cos(ws*t))/r3^3-ms*cos(ws*t)/rho^2;
    dxdt(4)=dxdt(4)-ms*(y-rho*sin(ws*t))/r3^3-ms*sin(ws*t)/rho^2;
    
    %Jacobians
    df1dx=1-(1-mu)/r1^3+3*(1-mu)*(x+mu)^2/r1^5-mu/r2^3+3*mu*(x+mu-1)^2/r2^5;
    df1dy=3*(1-mu)*(x+mu)*y/r1^5+3*mu*(x+mu-1)*y/r2^5;
    df2dy=1-(1-mu)/r1^3+3*(1-mu)*y^2/r1^5-mu/r2^3+3*mu*y^2/r2^5;

    df1dx=df1dx+ms*(2*(x-rho*cos(ws*t))^2-(y-rho*sin(ws*t))^2)/r3^5;
    df1dy=df1dy+3*ms*(x-rho*cos(ws*t))*(y-rho*sin(ws*t))/r3^5;
    df2dy=df2dy+ms*(2*(y-rho*sin(ws*t))^2-(x-rho*cos(ws*t))^2)/r3^5;

    A=[0, 0, 1, 0;            
       0, 0, 0, 1;
       df1dx, df1dy, 0, 2;
       df1dy, df2dy, -2, 0];

    dPHI = A*PHI;

    dxdt(5:end) = reshape(dPHI,16,1);
end

function [xx_ECI] = rot_to_ECI(xx_rot, mu, t) 
% DESCRIPTION
%   This function converts a state vector from a rotating synodic frame 
%   (e.g., Earth-Moon rotating frame) to the Earth-Centered Inertial (ECI) frame 
%   at a given dimensionless time. The conversion accounts for the rotation of 
%   the frame relative to the inertial frame.
%
% INPUT
%   xx_rot     [1x4]       state vector in the rotating frame, structured as 
%                          [x_position, y_position, x_velocity, y_velocity]             [adimensional]
%   mu         [1x1]       adimensional mass parameter                                  [-]
%   t          [1x1]       adimensional time                                            [-]
%
% OUTPUT
%   xx_ECI     [1x4]       state vector in the ECI frame, structured as 
%                          [x_position, y_position, x_velocity, y_velocity]             [adimensional]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB
xx_ECI = zeros(1,4);

xx_ECI(1) = (xx_rot(1) + mu)*cos(t) - xx_rot(2)*sin(t);
xx_ECI(2) = (xx_rot(1) + mu)*sin(t) + xx_rot(2)*cos(t);
xx_ECI(3) = (xx_rot(3) - xx_rot(2))*cos(t) - (xx_rot(4) + xx_rot(1) + mu)*sin(t);
xx_ECI(4) = (xx_rot(3) - xx_rot(2))*sin(t) + (xx_rot(4) + xx_rot(1) + mu)*cos(t);

end

function delta_v = deltavfun(xx_i)
% DESCRIPTION
%   This function computes the total Δv required for a two-impulse transfer 
%   in the Planar Bicircular Restricted Four-Body Problem (PBRFBP), using the 
%   initial and final states integrated over a given time interval. The function 
%   integrates the trajectory between the specified initial and final times and 
%   evaluates the velocity increments at departure and arrival.
%
% INPUT
%   xx_i      [1x6]       vector containing the initial state and transfer times:
%                         [x_position, y_position, x_velocity, y_velocity, t_initial, t_final]  [-]
%
% OUTPUT
%   delta_v   [1x1]       total Δv required for the transfer                                    [-]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

mu = 0.0121505842699404;

Re = 6378.1366;            
hi = 167;             
DU = 3.84405000*10^5;  
ri = (Re+hi)/DU;
Rm = 1737.4;           
hf = 100;             
rf=(Rm+hf)/DU;
ti = xx_i(5);
tf = xx_i(6);

options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[~, xx] = ode113(@(t,x) PBRFBP(t,x,mu), [ti tf], [xx_i(1);xx_i(2);xx_i(3);xx_i(4)], options);
xx_f = xx(end, :);

delta_v = (sqrt( (xx_i(3) -xx_i(2))^2 + (xx_i(4) + xx_i(1) + mu)^2) - sqrt((1-mu)/ri)) + abs(sqrt((xx_f(3) - xx_f(2))^2 + (xx_f(4) + xx_f(1) + mu - 1)^2) - sqrt(mu/rf)); 

end

function [c, ceq] = nonlcon(xx_i)
% DESCRIPTION
%   This function defines the nonlinear equality and inequality constraints 
%   for a two-impulse transfer in the Planar Bicircular Restricted Four-Body Problem (PBRFBP). 
%
% INPUT
%   xx_i      [1x6]       vector containing the initial state and transfer times:
%                         [x_position, y_position, x_velocity, y_velocity, t_initial, t_final]  [-]
%
% OUTPUT
%   c         []          vector of nonlinear inequality constraints (empty in this case)
%   ceq       [4x1]       vector of nonlinear equality constraints      [-]
%    
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

mu = 0.0121505842699404;

Re = 6378.1366;            
hi = 167;             
DU = 3.84405000*10^5;  
ri = (Re+hi)/DU;
Rm = 1737.4;           
hf = 100;             
rf=(Rm+hf)/DU;
ti = xx_i(5);
tf = xx_i(6);

ceq = zeros(4,1);

options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[~, xx] = ode113(@(t,x) PBRFBP(t,x,mu), [ti tf], [xx_i(1);xx_i(2);xx_i(3);xx_i(4)], options);
xx_f = xx(end, :);

ceq(1) = (xx_i(1) + mu)^2 + xx_i(2)^2 - ri^2;
ceq(2) = (xx_i(1)+mu)*(xx_i(3) - xx_i(2)) + xx_i(2)*(xx_i(4) + xx_i(1) + mu);
ceq(3) = (xx_f(1) + mu -1)^2 + xx_f(2)^2 - rf^2;
ceq(4) = (xx_f(1)+mu-1)*(xx_f(3) - xx_f(2)) + xx_f(2)*(xx_f(4) + xx_f(1) + mu-1);

c = [];
end

function h = circle(x,y,r,color) 
% DESCRIPTION
%   This function plots a circle with a specified center, radius, and optional color.
%
% INPUT
%   x       [1x1]       x-coordinate of the circle center                                   [-]
%   y       [1x1]       y-coordinate of the circle center                                   [-]
%   r       [1x1]       radius of the circle                                                [-]
%   color   [char]      (optional) color and line style specification for the plot command  [-]
%
% OUTPUT
%   h       [1x1]       command for plot                                                    [-]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

hold on
th = 0:0.01:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
if nargin>3
   h = plot(xunit, yunit,color);
else 
    h = plot(xunit, yunit);
end
end

function [dv, graddv] = deltavfun_grad(xx0)
% DESCRIPTION
%   This function computes the total Δv required for a PCR3BP
%   transfer between two specified positions at given times.
%   It also evaluates the gradient of the Δv with respect to the initial conditions 
%   and times, using state transition matrices (STM) for sensitivity analysis.
%
% INPUT
%   xx0     [6x1]   Vector of optimization variables:
%                   [xi; yi; vxi; vyi; ti; tf]
%                   where xi, yi, vxi, vyi are initial positions and velocities, and 
%                   ti, tf are initial and final times                                       [-]
%
% OUTPUT
%   dv      [1x1]   Total required Δv for the transfer                                       [-]
%   graddv  [6x1]   Gradient of the total Δv with respect to the decision variables          [-]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB
    % Constant input initialization
    mu = 0.0121505842699404;

    % Initialize initial point, times and target radius
    xi = xx0(1); yi = xx0(2);
    vxi = xx0(3); vyi = xx0(4);
    ti = xx0(5); tf = xx0(6);
    Re = 6378.1366;            
    hi = 167;             
    DU = 3.84405000*10^5;  
    ri = (Re+hi)/DU;
    Rm = 1737.4;           
    hf = 100;             
    rf=(Rm+hf)/DU;

    x0 = [xi; yi; vxi; vyi];
    graddv = zeros(6,1);

    % Initialize State Transition Matrix at t0 and append to initial 
    % conditions the conditions for the STM
    Phi0 = eye(4); x0Phi0 = [x0; Phi0(:)];

    % Propagate the initial condition
    options = odeset('reltol', 1e-12, 'abstol', 1e-12);
    [~, xxf] = ode113(@(t,x) PBRFBP_STM(t,x,mu), [ti, tf], x0Phi0, options);

    % Initialize final point
    xf = xxf(end,1); yf = xxf(end,2);
    vxf = xxf(end,3); vyf = xxf(end,4);

    % Reshape STM
    PHIf = reshape(xxf(end,5:end),4,4);

    % Delta velocity evaluation
    vi = sqrt((vxi - yi)^2 + (vyi + xi + mu)^2);
    vf = sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2);
    dvi = vi - sqrt((1-mu)/ri);
    dvf = vf - sqrt(mu/rf);

    % Objective function and gradient definition
    dv = dvi + dvf;
    if nargout > 1
        graddv(1:4) = 1 / vi * [vyi+xi+mu; yi-vxi; -yi+vxi; vyi+xi+mu] + PHIf' * 1/ vf * [vyf+xf+mu-1; yf-vxf; -yf+vxf; vyf+xf+mu-1];
        graddv(5) = (1/ vf * [vyf+xf+mu-1; yf-vxf; -yf+vxf; vyf+xf+mu-1])' * (-PHIf)* PBRFBP(ti,xx0(1:4),mu);
        graddv(6) = (1/ vf * [vyf+xf+mu-1; yf-vxf; -yf+vxf; vyf+xf+mu-1])' * PBRFBP(tf,[xf;yf;vxf;vyf],mu);
    end
end

function [c, ceq, gradc, gradceq] = nonlcon_grad(xx) 
% DESCRIPTION
%   This function computes the nonlinear equality constraints and their gradients 
%   for a PCR3BP transfer. 
%
% INPUT
%   xx      [6x1]   Vector of decision variables:
%                   [xi; yi; vxi; vyi; ti; tf]
%                   where xi, yi, vxi, vyi are initial positions and velocities, and 
%                   ti, tf are initial and final times                                            [-]
%
% OUTPUT
%   c           []      Inequality constraints (empty)                                            [-]
%   ceq         [4x1]   Nonlinear equality constraints                                            [-]
%   gradc       []      Gradient of inequality constraints (empty)                                [-]
%   gradceq     [6x4]   Gradient of equality constraints with respect to optimization variables   [-]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

    % Constant input initialization
    mu = 0.0121505842699404;

    % Initialize initial point, times and target radius
    xi = xx(1); yi = xx(2);
    vxi = xx(3); vyi = xx(4);
    ti = xx(5); tf = xx(6);
    Re = 6378.1366;            
    hi = 167;             
    DU = 3.84405000*10^5;  
    ri = (Re+hi)/DU;
    Rm = 1737.4;           
    hf = 100;             
    rf=(Rm+hf)/DU;

    x0 = [xi; yi; vxi; vyi];

    % Initialize State Transition Matrix at t0 and append to initial 
    % conditions the conditions for the STM
    Phi0 = eye(4); x0Phi0 = [x0; Phi0(:)];

    % Propagate the initial condition
    options = odeset('reltol', 1e-12, 'abstol', 1e-12);
    [~, xxf] = ode78(@(t,x) PBRFBP_STM(t,x,mu), [ti, tf], x0Phi0, options);
    
    % Initialize final point
    xf = xxf(end,1); yf = xxf(end,2);
    vxf = xxf(end,3); vyf = xxf(end,4);

    % Reshape STM
    PHIf = reshape(xxf(end,5:end),4,4);

    % Non-Linear Constrains definition
    c = [] ;
    ceq = [(xi + mu)^2 + yi^2 - ri^2, (xi + mu)*(vxi - yi) + yi*(vyi + xi + mu), (xf + mu - 1)^2 + yf^2 - rf^2, (xf + mu - 1)*(vxf - yf) + yf*(vyf + xf + mu - 1)];

    if nargout > 2
        % Non-Linear Constrains Gradient definition
        gradc = [];
        gradceq = [2*(xi+mu),2*yi,0,0,0,0;
            vxi,vyi,xi+mu,yi,0,0;
            [2*mu + 2*xf - 2, 2*yf, 0, 0]*PHIf, -[2*mu + 2*xf - 2, 2*yf, 0, 0]*PHIf*PBRFBP(ti, xx(1:4),mu), [2*mu + 2*xf - 2, 2*yf, 0, 0]*PBRFBP(tf,[xf;yf;vxf;vyf],mu);
            [vxf, vyf, mu + xf - 1, yf]*PHIf, -[vxf, vyf, mu + xf - 1, yf]*PHIf*PBRFBP(ti, xx(1:4),mu), [vxf, vyf, mu + xf - 1, yf]*PBRFBP(tf,[xf;yf;vxf;vyf],mu)];
        gradceq = gradceq.' ;
    end
end

function [dv, graddv] = deltavfun_MS_grad(xx0)
% DESCRIPTION
%   This function computes the total required Δv for a multiple shooting
%   transfer with 4 nodes in the PBRFBP.
%   based on initial and final state vectors. It also computes the gradient of 
%   the Δv with respect to the optimization variables.
% INPUT
%   xx0    [18x1]   State vector containing initial and final position and velocity 
%                    components for N shots, adimensionalized                      [-]
% OUTPUT
%   dv     [1x1]     Total required delta-v for the trajectory, adimensionalized   [-]
%   graddv [18x1]   Gradient of delta-v with respect to the state vector           [-]
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

    % Constant input initialization
    mu = 0.0121505842699404;    
    % Number of shots
    N = 4;
    % Initialize initial points, constants and target radius
    xi = xx0(1); yi = xx0(2); vxi = xx0(3); vyi = xx0(4); 
    xf = xx0(4*N-3); yf = xx0(4*N-2); vxf = xx0(4*N-1); vyf = xx0(4*N);
    Re = 6378.1366;            
    hi = 167;             
    DU = 3.84405000*10^5;  
    ri = (Re+hi)/DU;
    Rm = 1737.4;           
    hf = 100;             
    rf=(Rm+hf)/DU;

    % Delta velocity evaluation
    vi = sqrt((vxi - yi)^2 + (vyi + xi + mu)^2);
    vf = sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2);
    dvi = vi - sqrt((1-mu)/ri);
    dvf = vf - sqrt(mu/rf);

    % Objective function definition
    dv = dvi + dvf;

    % Gradient definition
    graddv = zeros(4*N+2,1);
    if nargout > 1
        graddv(1:4) = 1 / vi * [vyi+xi+mu; yi-vxi; -yi+vxi; vyi+xi+mu];
        graddv((4*N-3):4*N) = 1/ vf * [vyf+xf+mu-1; yf-vxf; -yf+vxf; vyf+xf+mu-1];
    end

end

function [c, ceq, gradc, gradceq] = nonlcon_MS_grad(xx0)
% DESCRIPTION
%   This function evaluates the nonlinear inequality and equality constraints, 
%   along with their gradients, for a 4 nodes multiple shooting trajectory optimization problem 
%   in the PBRFBP.
% INPUT
%   xx0     [18x1]   State vector containing the concatenated position, velocity, 
%                         and initial/final times for N shots, adimensionalized      [-]
% OUTPUT
%   c       [8x1]      Vector of nonlinear inequality constraints                    [-]
%   ceq     [16x1]      Vector of nonlinear equality constraints                     [-]
%   gradc   [18x8]  Gradient of inequality constraints with respect to xx0           [-]
%   gradceq [18x16]  Gradient of equality constraints with respect to xx0            [-]
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

    % Constant input initialization
    mu = 0.0121505842699404;
    % Number of shots
    N = 4;
    % Initialize initial point, times and target radius
    xx = zeros(4,N);
    for i = 1 : N
        xx(:,i) = xx0(4*i-3 : 4*i);
    end
    ti = xx0(end-1); tf = xx0(end);
    Re = 6378.1366;            
    hi = 167;             
    DU = 3.84405000*10^5;  
    ri = (Re+hi)/DU;
    Rm = 1737.4;           
    hf = 100;             
    rf=(Rm+hf)/DU;
    Rm = 1737.4/DU;
    Re = 6378.1366/DU; 
    % Initialize State Transition Matrix at t0, vector of solutions and options for integration
    PHI0 = eye(4); xxf = zeros(4,N-1); PHIf = zeros(4,4,N-1);
    options = odeset('reltol', 1e-12, 'abstol', 1e-12);
    t_vect = linspace(ti,tf,N);
    for i = 1 : N-1
        % Propagate the initial condition
        [~, xf] = ode78(@(t,x) PBRFBP_STM(t,x,mu), [t_vect(i), t_vect(i+1)],[xx(:,i); PHI0(:)], options);
        % Saving final integred points
        xxf(:,i) = xf(end,1:4)';
        % Reshape STM
        PHIf(:,:,i) = reshape(xf(end,5:end),4,4);
    end

    % Non-Linear Constrains definition
    c = zeros(2*N,1);
    for i = 1 : N
        c((2*i-1):2*i) = [Re^2 - (xx(1,i) + mu)^2 - xx(2,i)^2 ;
                            Rm^2 - (xx(1,i)+mu-1)^2 - xx(2,i)^2];
    end
    % c(end) = ti-tf ;
    ceq = zeros(4*N,1);
    for i = 1:N-1
        ceq((4*i-3):4*i) = xxf(:,i) - xx(:,i+1);
    end
    ceq(end-3:end) = [(xx(1,1)+mu)^2 + xx(2,1)^2 - ri^2;
                    (xx(1,1)+mu)*(xx(3,1)-xx(2,1)) + xx(2,1)*(xx(4,1) + xx(1,1) + mu);
                    (xx(1,N)+mu-1)^2 + xx(2,N)^2 - rf^2;
                    (xx(1,N)+mu-1)*(xx(3,N)-xx(2,N)) + xx(2,N)*(xx(4,N) + xx(1,N) + mu - 1)];


    if nargout > 2
        
        % Non-Linear Constrains Gradient definition
        const = mu;
        gradc = zeros(2*N,4*N+2);
        for i = 1 : N
            gradc((2*i-1):2*i,(4*i-3):4*i) = [-2*(xx(1,i)+mu), -2*xx(2,i), 0, 0;
                                                -2*(xx(1,i)+mu-1), -2*xx(2,i), 0, 0];
        end
        gradceq = zeros(4*N,4*N+2); 
        for i = 1 : N-1
            gradceq(4*i-3:4*i,4*i-3:4*i+4) = [PHIf(:,:,i), -eye(4)]; 
        end
        gradceq(1:4,17) = -PHIf(:,:,1)*PBRFBP(t_vect(1),xx(:,1),mu) + PBRFBP(t_vect(2),xxf(:,1),mu)*(N-2)/(N-1);
        gradceq(1:4,end) =  PBRFBP(t_vect(2),xxf(:,1),mu)/(N-1);
        gradceq(5:8,17) = -PHIf(:,:,2)*PBRFBP(t_vect(2),xx(:,2),mu)*(N-2)/(N-1) + PBRFBP(t_vect(3),xxf(:,2),mu)*(N-3)/(N-1);
        gradceq(5:8,end) =  -PHIf(:,:,2)*PBRFBP(t_vect(2),xx(:,2),mu)/(N-1) + 2*PBRFBP(t_vect(3),xxf(:,2),mu)/(N-1);
        gradceq(9:12,17) = -PHIf(:,:,3)*PBRFBP(t_vect(3),xx(:,3),const)*(N-3)/(N-1);
        gradceq(9:12,end) = -2*PHIf(:,:,3)*PBRFBP(t_vect(3),xx(:,3),const)/(N-1) + 3*PBRFBP(t_vect(4),xxf(:,3),const)/(N-1);



        gradceq(end-3:end-2,1:4) = [2*(xx(1,1)+mu), 2*xx(2,1), 0, 0; %R1
                                    xx(3,1), xx(4,1), xx(1,1)+mu, xx(2,1)];
        gradceq(end-1:end,end-5:end-2) = [2*(xx(1,N)+mu-1), 2*xx(2,N), 0, 0; %R2
                                            xx(3,N), xx(4,N), xx(1,N)+mu-1, xx(2,N)];
         gradc = gradc';
         gradceq = gradceq'; 
    end
end

function [dxdt] = nbody(tt, xx, bodies, frame, center)
% DESCRIPTION
%   This function evaluates the right-hand side of a N-body propagator 
%   for orbital dynamics, computing the time derivative of a spacecraft’s 
%   Cartesian state vector under the influence of gravitational 
%   accelerations from multiple celestial bodies. The integration can be performed 
%   in either the 'J2000' or 'ECLIPJ2000' reference frames, with a specified 
%   integration center.
%
% INPUT
%   tt      [1x1]       ephemeris time (ET SPICE), seconds past J2000 (TDB)             [s]
%   xx      [6x1]       Cartesian state vector of the object with respect to the 
%                       desired integration center, structured as [position(3); velocity(3)]  [km] [km/s]
%   bodies  [1xn]       cell array of strings containing the names of celestial bodies 
%                       involved in the propagation (as SPICE-recognized names)         [-]
%   frame   [char]      string specifying the reference frame, either 'J2000' or 
%                       'ECLIPJ2000'                                                    [-]
%   center  [char]      string specifying the name of the central body in the propagation 
%                       (must be contained within 'bodies')                             [-]
%
% OUTPUT
%   dxdt    [6x1]       time derivative of the state vector containing the velocity 
%                       and acceleration of the object                                  [km/s] [km/s^2]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB
    
    if not( strcmpi(frame, 'ECLIPJ2000') || strcmpi(frame, 'J2000') )
        msg = 'Invalid integration reference frame, select either J2000 or ECLIPJ2000';
        error(msg);
    end

    if not(any(strcmp(center,bodies)))
        msg = 'Invalid center, select a bodies array containing the center body';
        error(msg);
    end
    
    % Gravitational Parameter initialization
    GM = zeros(length(bodies),1);
    for i = 1 : length(bodies)
        if strcmp(center,bodies(i))
            GM0 = cspice_bodvrd(bodies(i),'GM',1);
        end
        GM(i) = cspice_bodvrd(bodies(i),'GM',1);
    end

    % Initialize right-hand-side
    dxdt = zeros(6,1);

    % Position detivative is object's velocity
    dxdt (1:3) = xx(4:6);

    % Compute contribution to acceleration of central body
    dxdt(4:6) = - GM0 * xx(1:3) / norm(xx(1:3))^3 ;

    % Loop over all bodies (except GM0)
    for i = 1:length(GM)
        if not(strcmp(center,bodies(i)))

            % Retrieve position and velocity of i-th celestial body wrt desired
            % center in non-inertial frame:
            rv_b0_body = cspice_spkezr(bodies(i), tt, frame, 'NONE', center);

            % Extract object position wrt. i-th celestial body
            rr_body_obj = -rv_b0_body(1:3) + xx(1:3) ;

            % Compute non-inertial terms
            rho = rv_b0_body(1:3);
            d = rr_body_obj ;
            q = dot(xx(1:3),xx(1:3)-2*rho)/dot(rho,rho);
            f = q * (3 + 3*q + q^2) / (1+ (1 + q)^1.5);

            % Compute their contribution to acceleration:
            aa_grav = - GM(i) * 1/norm(d)^3 *(xx(1:3) + rho*f);

            % Sum up acceleration to right-hand-side
            dxdt(4:6) = dxdt(4:6) + aa_grav ;
        end
    end
end
