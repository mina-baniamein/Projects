%Spacecraft Guidance and Navigation
%Assignment #1
%Author: Mina Baniamein

%% Point 1
clear; close all; clc;
% Inital guesses for collinear points chosen graphically from dUdx plot 

mu = 0.012150;
earth_x = -mu; % x coordinate of Earth 
moon_x = 1-mu; % x coordinate of Moon
dUdx = @(x) x-(1-mu).*(x+mu)./(abs(x+mu)).^3-mu.*(x+mu-1)./(abs(x+mu-1)).^3;

figure
hold on
plot(earth_x, 0, 'ko', 'MarkerSize', 10,'MarkerFaceColor', '#0019D1'); % Earth
plot(moon_x, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', '#7C8886'); % Moon
x1 = -2:0.01:-0.21 ; x2 = 0.1:0.01:0.97 ; x3 = 1.01:0.01:2; 
plot(x1, dUdx(x1),'b', x2, dUdx(x2),'b', x3, dUdx(x3),'b') ; axis([-2 2 -30 30]) ;
grid on
xlabel('x [-]');
ylabel('$\mathbf{\frac{\partial U}{\partial x}}$','Interpreter','Latex');
title('$\mathbf{\frac{\partial U}{\partial x}}$ plot','Interpreter','Latex');
ylim([-20, 20]); 
legend('Earth','Moon', '$\frac{\partial U}{\partial x}$','location', 'northwest','interpreter','latex')
title('Plot for collinear points guess')

%L1 x coordinate computation
options = optimoptions('fsolve','OptimalityTolerance',1e-12);
x_L1 = fsolve(dUdx, -0.94, options);
C_L1 = jacobi_const(x_L1,0,0,0,0,0,mu);
disp(x_L1); 

%L2 x coordinate computation
x_L2 = fsolve(dUdx, 0.8, options);
C_L2 = jacobi_const(x_L2,0,0,0,0,0,mu);

%L3 x coordinate computation
x_L3 = fsolve(dUdx, 1.2, options);
C_L3 = jacobi_const(x_L3,0,0,0,0,0,mu);

% Inital guesses for triangular points chosen graphically from cirlces 
% intersections in plot 

figure
hold on
grid on
axis equal
plot(earth_x, 0, 'ko', 'MarkerSize', 10,'MarkerFaceColor', '#0019D1'); 
plot(moon_x, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', '#7C8886'); 
circle(earth_x,0,1)
circle(moon_x,0,1)
xlabel('x [-]');
ylabel('y [-]');
legend('Earth','Moon')
title('Plot for triangular points guess')

F = @(x) [(x(1) + mu).^2 + x(2).^2-1;
          (x(1) -1 + mu).^2 + x(2).^2-1];
% L4 coordinates computation
x0_L4 = [0.5;0.86];

[r_L4, ~] = fsolve(F,x0_L4,options);
C_L4 = jacobi_const(r_L4(1),r_L4(2), 0, 0, 0, 0, mu);

% L5 coordinates computation
x0_L5 = [0.5;-0.86];

[r_L5, ~] = fsolve(F,x0_L5,options);
C_L5 = jacobi_const(r_L5(1),r_L5(2), 0, 0, 0, 0, mu);

%Results
clc
disp('Coordinates of Earth-Moon lagrangian points on xy plane:');
disp(['L1: (', num2str(x_L1, '%.10f'),',0)']);
disp(['L2: (', num2str(x_L2, '%.10f'),',0)']);
disp(['L3: (', num2str(x_L3, '%.10f'),',0)']);
disp(['L4: (', num2str(r_L4(1), '%.10f'),',', num2str(r_L4(2), '%.10f'),')']);
disp(['L5: (', num2str(r_L5(1), '%.10f'),',', num2str(r_L5(2), '%.10f'),')']);

disp('Jacobi constant value of Earth-Moon lagrangian points:');
disp(['C L1: ', num2str(C_L1, '%.10f')]);
disp(['C L2: ', num2str(C_L2, '%.10f')]);
disp(['C L3: ', num2str(C_L3, '%.10f')]);
disp(['C L4: ', num2str(C_L4, '%.10f')]);
disp(['C L5: ', num2str(C_L5, '%.10f')]);

figure
hold on
grid on
axis equal
plot(earth_x, 0, 'ko', 'MarkerSize', 10,'MarkerFaceColor', '#0019D1'); 
plot(moon_x, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', '#7C8886'); 
plot(x_L1,0,'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
plot(x_L2,0,'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
plot(x_L3,0,'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'c');
plot(r_L4(1),r_L4(2),'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'm');
plot(r_L5(1),r_L5(2),'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
plot(x1, dUdx(x1),'b', x2, dUdx(x2),'b', x3, dUdx(x3),'b') ; axis([-2 2 -30 30]) ;
ylim([-2, 2]);
circle(earth_x,0,1);
circle(moon_x,0,1);
legend('Earth','Moon','L1','L2','L3','L4','L5','Interpreter','Latex','location', 'northwest');
xlabel('x [-]');
ylabel('y [-]');
title('Lagrangian points location')

%% Point 2
clc; close all;

mu = 0.012150;

% Initial conditions
x0  = 1.068792441776;
y0  = 0;
z0  = 0.071093328515;
vx0 = 0;
vy0 = 0.319422926485;
vz0 = 0;

phi0 = reshape(eye(6),36,1);

xx0 = [x0;y0;z0;vx0;vy0;vz0;phi0];

C_0 = jacobi_const(x0,y0,z0,vx0,vy0,vz0, mu);

%%% Correction of initial states and correct orbit propagation:


% Define solver's parameters
Nmax    = 50;     % Maximum number of iterations
tol     = 1e-10;    % Desired tolerance

% Find correct initial states and propagate the orbit
opt = odeset('RelTol',2.22045e-14,'AbsTol',1e-20,'Events',@(x,y) x_axis_crossing(x,y,1));

iter = 0;           % Iteration counter
errors = ones(4,1); % Initial guess of error on vx an vz (higher than tolerance)


% start Newton Method
while max(abs(errors)) > tol && iter < Nmax

    %%% Compute the new initial conditions

    if iter % if this is not the first iteration, compute the new guesses

        % Select the only elements of the STM useful to correct the previous
        % guess
        phi = reshape(final_states(7:end),6,6);
        x = final_states(1);
        y = final_states(2);
        z = final_states(3);
        r1 = sqrt((x + mu)^2 + y^2+z^2);
        r2 = sqrt((x + mu - 1)^2 + y^2+z^2);
        dUdx = x - (1-mu)/r1^3*(mu+x) + mu/r2^3*(1-mu-x);
        dUdz = - (1-mu)/r1^3*z -mu/r2^3*z;

        phi_p = [phi(2,1),phi(2,3),phi(2,5),final_states(5); 
                 phi(4,1), phi(4,3), phi(4,5), 2*final_states(5) + dUdx; 
                 phi(6,1), phi(6,3), phi(6,5), dUdz; 
                 2*x + ((2*mu + 2*x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) - (mu*(2*mu + 2*x - 2))/((mu + x - 1)^2 + y^2 + z^2)^(3/2), (2*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) - (2*mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2), -2*final_states(5), 0;];
        
        % Compute the new guess
        delta_x0 = phi_p\[final_states(2);final_states(4);final_states(6); 3.09 - C_0];
        xx0([1;3;5]) = xx0([1;3;5]) - delta_x0([1;2;3]); %added damping
        % Compute Jacobian constant
        C_0 = jacobi_const(xx0(1), xx0(2), xx0(3), xx0(4), xx0(5), xx0(6), mu);
        

    end

    % Propagate the orbit with the new initial conditions
    [~,states] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,100],xx0,opt);
    final_states = states(end,:);

    % Compute the error between the reference vx_final and vz_final (=0)
    % and the ones obtained propagating the dynamics
    errors = [final_states(2); final_states(4); final_states(6); 3.09 - C_0];

    %increase iteration counter
    iter = iter + 1;
    
end

sol = xx0;
rv_sol = sol(1:6);
errors;
C_sol = jacobi_const(xx0(1), xx0(2), xx0(3), xx0(4), xx0(5), xx0(6), mu);

disp('Corrected initial position and velocity:');
disp(['x_0:  ', num2str(rv_sol(1), '%.10f')]);
disp(['y_0:  ', num2str(rv_sol(2), '%.10f')]);
disp(['z_0:  ', num2str(rv_sol(3), '%.10f')]);
disp(['vx_0: ', num2str(rv_sol(4), '%.10f')]);
disp(['vy_0: ', num2str(rv_sol(5), '%.10f')]);
disp(['vz_0: ', num2str(rv_sol(6), '%.10f')]);
disp(' ');
disp(['Number of iterations: ', num2str(iter)]);

% Halo orbit plot
[~,states] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,100],sol,opt);
[~,back_states] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,-100],sol,opt);

figure
plot3(states(:,1),states(:,2),states(:,3))
hold on
grid on
plot3(back_states(:,1),back_states(:,2),back_states(:,3))
plot(moon_x, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', '#7C8886'); 
plot(x_L3,0,'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'c');
legend('Halo orbit - forward prop.','Halo orbit - backward prop.','Moon','L3','location', 'northwest');
xlabel('x [-]');
ylabel('y [-]');
zlabel('z [-]')
title('C = 3.09 halo orbit')


%% Point 3
clc; close all;
C_f = 3.04;
C_i = 3.09;
C_halos = linspace(C_i,C_f,6); %constants list of halo family
halos = zeros(6,6);  %storing vector for initial conditions of halo orbits propagations

xx0 = sol;
C_0 = jacobi_const(xx0(1), xx0(2), xx0(3), xx0(4), xx0(5), xx0(6), mu);

for i = 1:length(C_halos)
    iter = 0;  % Reset iteration counter
    errors = ones(4,1);  % Reset error for each C_halos(i)

% start Newton Method
while max(abs(errors)) > tol && iter < Nmax
    %%% Compute the new initial conditions

    if iter % if this is not the first iteration, compute the new guesses

        % Select the only elements of the STM useful to correct the previous
        % guess
        phi = reshape(final_states(7:end),6,6);
        x = final_states(1);
        y = final_states(2);
        z = final_states(3);
        r1 = sqrt((x + mu)^2 + y^2+z^2);
        r2 = sqrt((x + mu - 1)^2 + y^2+z^2);
        dUdx = x - (1-mu)/r1^3*(mu+x) + mu/r2^3*(1-mu-x);
        dUdz = - (1-mu)/r1^3*z -mu/r2^3*z;

        phi_p = [phi(2,1),phi(2,3),phi(2,5),final_states(5); 
                 phi(4,1), phi(4,3), phi(4,5), 2*final_states(5) + dUdx; 
                 phi(6,1), phi(6,3), phi(6,5), dUdz; 
                 2*x + ((2*mu + 2*x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) - (mu*(2*mu + 2*x - 2))/((mu + x - 1)^2 + y^2 + z^2)^(3/2), (2*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) - (2*mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2), -2*final_states(5), 0;];
        
        % Compute the new guess
        delta_x0 = phi_p\[final_states(2);final_states(4);final_states(6); C_halos(i) - C_0];
        xx0([1;3;5]) = xx0([1;3;5]) - delta_x0([1;2;3]); 
        C_0 = jacobi_const(xx0(1), xx0(2), xx0(3), xx0(4), xx0(5), xx0(6), mu);
        %xx0([1;3;5]) = guess;

    end

    % Propagate the orbit with the new initial conditions
    [~,states] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,100],xx0,opt);
    final_states = states(end,:);

    % Compute the error between the reference vx_final and vz_final (=0)
    % and the ones obtained propagating the dynamics
    errors = [final_states(2); final_states(4); final_states(6); C_halos(i) - C_0];

    %increase iteration counter
    iter = iter + 1;
    
end
xx0_halo = xx0;
halos(i,:) = xx0_halo(1:6)';
end

halos;
% Verification of numerical continuation being respected
C2 = jacobi_const(halos(2,1), halos(2,2), halos(2,3), halos(2,4), halos(2,5), halos(2,6), mu);
C3 = jacobi_const(halos(3,1), halos(3,2), halos(3,3), halos(3,4), halos(3,5), halos(3,6), mu);
C4 = jacobi_const(halos(4,1), halos(4,2), halos(4,3), halos(4,4), halos(4,5), halos(4,6), mu);
C5 = jacobi_const(halos(5,1), halos(5,2), halos(5,3), halos(5,4), halos(5,5), halos(5,6), mu);
C6 = jacobi_const(halos(6,1), halos(6,2), halos(6,3), halos(6,4), halos(6,5), halos(6,6), mu);

%%
disp('Resulting Jacobi constants from computed initial conditions:');
disp(['C2: ', num2str(C2)]);
disp(['C3: ', num2str(C3)]);
disp(['C4: ', num2str(C4)]);
disp(['C5: ', num2str(C5)]);
disp(['C6: ', num2str(C6)]);
disp(' ');
disp('Corrected initial position and velocity with C = 3.04:');
disp(['x_0:  ', num2str(halos(end,1), '%.10f')]);
disp(['y_0:  ', num2str(halos(end,2), '%.10f')]);
disp(['z_0:  ', num2str(halos(end,3), '%.10f')]);
disp(['vx_0: ', num2str(halos(end,4), '%.10f')]);
disp(['vy_0: ', num2str(halos(end,5), '%.10f')]);
disp(['vz_0: ', num2str(halos(end,6), '%.10f')]);
disp(' ');
disp(['Number of iterations: ', num2str(iter)]);
%% 

% Family of Halo orbits plot
[~,states_C1] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,100],[halos(1,:)'; phi0],opt);
[~,back_states_C1] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,-100],[halos(1,:)'; phi0],opt);
[~,states_C2] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,100],[halos(2,:)'; phi0],opt);
[~,back_states_C2] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,-100],[halos(2,:)'; phi0],opt);
[~,states_C3] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,100],[halos(3,:)'; phi0],opt);
[~,back_states_C3] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,-100],[halos(3,:)'; phi0],opt);
[~,states_C4] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,100],[halos(4,:)'; phi0],opt);
[~,back_states_C4] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,-100],[halos(4,:)'; phi0],opt);
[~,states_C5] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,100],[halos(5,:)'; phi0],opt);
[~,back_states_C5] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,-100],[halos(5,:)'; phi0],opt);
[~,states_C6] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,100],[halos(6,:)'; phi0],opt);
[~,back_states_C6] = ode113(@(t,xx) CR3BP_STM_3D(t,xx,mu),[0,-100],[halos(6,:)'; phi0],opt);

figure
plot3(states_C1(:,1),states_C1(:,2),states_C1(:,3),'b')
hold on
grid on
plot3(back_states_C1(:,1),back_states_C1(:,2),back_states_C1(:,3),'b')

plot3(states_C2(:,1),states_C2(:,2),states_C2(:,3),'r')
hold on
grid on
plot3(back_states_C2(:,1),back_states_C2(:,2),back_states_C2(:,3),'r')

plot3(states_C3(:,1),states_C3(:,2),states_C3(:,3),'y')
hold on
grid on
plot3(back_states_C3(:,1),back_states_C3(:,2),back_states_C3(:,3),'y')

plot3(states_C4(:,1),states_C4(:,2),states_C4(:,3),'g')
hold on
grid on
plot3(back_states_C4(:,1),back_states_C4(:,2),back_states_C4(:,3),'g')

plot3(states_C5(:,1),states_C5(:,2),states_C5(:,3),'c')
hold on
grid on
plot3(back_states_C5(:,1),back_states_C5(:,2),back_states_C5(:,3),'c')

plot3(states_C6(:,1),states_C6(:,2),states_C6(:,3),'m')
hold on
grid on
plot3(back_states_C6(:,1),back_states_C6(:,2),back_states_C6(:,3),'m')

plot(moon_x, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', '#7C8886'); 
plot(x_L3,0,'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'c');
legend('','', '','','','', '','', '','','','', 'Moon', 'L3', 'location', 'northwest');
xlabel('x [-]');
ylabel('y [-]');
zlabel('z [-]')
title('From C = 3.09 to C = 3.04 halo orbit family')







%% functions

function h = circle(x,y,r) 
% DESCRIPTION
%   This function plots a circle of desired center and radius 
% INPUT
%   x    x coordinate of the circle center  [-]
%   y    y coordinate of the circle center  [-]
%   r    Radius of the circle               [-]
% OUTPUT
%   h: Plot command
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

hold on
th = 0:0.01:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end

function C = jacobi_const(x,y,z,vx,vy,vz, mu) %calculates Jacobian constant 
% DESCRIPTION
%   This function computes the Jacobi constant
% INPUT
%   x  [1x1]     x coordinate of the adimensional position vector [-]
%   y  [1x1]     y coordinate of the adimensional position vector [-]
%   z  [1x1]     z coordinate of the adimensional position vector [-]
%   vx [1x1]     x component of the adimensional velocity vector  [-]
%   vy [1x1]     y component of the adimensional velocity vector  [-]
%   vz [1x1]     z component of the adimensional velocity vector  [-]
%   mu [1x1]     adimensionalized mass parameter                  [-]
% OUTPUT
%   h   Plot command
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

r1 = sqrt((x+mu).^2 + y.^2 + z.^2);
r2 = sqrt((x + mu -1).^2 + y.^2 + z.^2);
v_squared = vx.^2 + vy.^2 + vz.^2;
C = (0.5*(x.^2 + y.^2) + (1-mu)./r1 + mu./r2 +0.5.*mu.*(1-mu)).*2 - v_squared;
end

function dx = CR3BP_STM_3D(~,xx,mu)
% DESCRIPTION
%   This function contains the RHS for solving the CR3BP differential
%   equations and computing the state transition matrix
% INPUT
%   xx  [6x1]     state vector containing adimensionalized position and velocity vectors [-]
%   mu  [1x1]     adimensionalized mass parameter                                        [-]
% OUTPUT
%   dx  [42x1]    vector containing derivative of the state vector on the                [-] 
%                 first 6 rows, and elements of the 6x6 STM on the
%                 other 36 rows
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

% s/c states
x = xx(1);
y = xx(2);
z = xx(3);
vx = xx(4);
vy = xx(5);
vz = xx(6);

% Compute distances from bodies 1 and 2
r1 = sqrt((x + mu)^2 + y^2+z^2);
r2 = sqrt((x + mu - 1)^2 + y^2+z^2);
% compute partial derivatives of scalar potential function U
dUdx = x - (1-mu)/r1^3*(mu+x) + mu/r2^3*(1-mu-x);
dUdy = y - (1-mu)/r1^3*y - mu/r2^3*y;
dUdz = - (1-mu)/r1^3*z -mu/r2^3*z;

% assemble right-hand side
dx = zeros(42,1);
dx(1:3) = [vx;vy;vz];
dx(4) = 2*vy + dUdx;
dx(5) = -2*vx + dUdy;
dx(6) = dUdz;


%%% Integration of the STM

% Put the STM PHI in matrix form
phi = reshape(xx(7:end),6,6);

% Jacobian of f
A = zeros(6);
A(1,4) = 1;
A(2,5) = 1;
A(3,6) = 1;
A(4,1) = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu + 2*x)*(mu + x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)) + 1;
A(4,2) = (3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
A(4,3) = (3*mu*z*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
A(4,5) = 2;
A(5,1) = (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));
A(5,2) = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1;
A(5,3) = (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
A(5,4) = -2;
A(6,1) = (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));
A(6,2) = (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
A(6,3) = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2);


%Compute the derivative of the STM
dphi = A*phi;

dx(7:end) = reshape(dphi,36,1);

end

function [value, isterminal, direction] = x_axis_crossing(~,xx,isTerminal)
% DESCRIPTION
%   Event function that arrest numerical solvers when the S/C crosses the x axis 
% INPUT
%   xx          [6x1]     state vector containing adimensionalized position and velocity vectors [-]
%   isTerminal  [1x1]     flag stopping event if set to 1                                        [-]
% OUTPUT
%   value       [1x1]     checks when the y coordinate value is 0                                [-] 
%   isterminal  [1x1]     set to input isTerminal to decide if arrest the numerical solver       [-]        
%   direction   [1x1]     set to zero in order to ignore crossing verse of the x axis            [-]             
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB
    value = xx(2);
    isterminal = isTerminal;
    direction = 0;
end
