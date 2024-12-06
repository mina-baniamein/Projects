%Spacecraft Guidance and Navigation
%Assignment #1
%Author: Mina Baniamein

%% Point 1
clear ; close all; clc;
%% 

mu = 0.012150;
rE = -mu; %Earth %potrebbero non essere necessarie queste variabili
rM = 1-mu; %Moon
dUdx = @(x) x-(1-mu).*(x+mu)./(abs(x+mu)).^3-mu.*(x+mu-1)./(abs(x+mu-1)).^3;
x_axis = linspace(-2,2, 100);
plot(x_axis, dUdx(x_axis))
grid on
ylim([-0.4 0.8])


%% 
%mettere guess
%L1
options = optimoptions('fsolve','OptimalityTolerance',1e-10,'Display','iter-detailed');
x_L1 = fsolve(dUdx, -0.94, options)
%C_L1 = jacobiLAG(x_L1,0,mu)
 

%L2
x_L2 = fsolve(dUdx, 0.8, options)
%C_L2 = jacobiLAG(x_L2,0,mu)

%L3
x_L3 = fsolve(dUdx, 1.2, options)
%C_L3 = jacobiLAG(x_L3,0,mu)

%% 

%Earth and Moon position and circles intersections 

%plot for initial guess
figure
hold on
circle(-mu,0,1)
circle(1-mu,0,1)

F = @(x) [(x(1) + mu).^2 + x(2).^2-1;
          (x(1) -1 + mu).^2 + x(2).^2-1];

x0_L4 = [0.5;0.86];

[r_L4, toll] = fsolve(F,x0_L4,options)
C_L4 = jacobi_const(r_L4(1),r_L4(2), 0, 0, 0, 0, mu)


x0_L5 = [0.5;-0.86];

[r_L5, toll] = fsolve(F,x0_L5,options)
C_L5 = jacobi_const(r_L5(1),r_L5(2), 0, 0, 0, 0, mu)
%% Point 2
clearvars; close all; clc;

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
        phi = reshape(final_states(7:end),6,6)
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
        C_0 = jacobi_const(xx0(1), xx0(2), xx0(3), xx0(4), xx0(5), xx0(6), mu);
        %xx0([1;3;5]) = guess;

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
xx0
iter
errors
%% 
jacobi_const(xx0(1), xx0(2), xx0(3), xx0(4), xx0(5), xx0(6), mu)







%% functions
%circle
function h = circle(x,y,r) %for circles graphing
hold on
th = 0:0.01:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end
function C = jacobi_const(x,y,z,vx,vy,vz, mu) %calculates Jacobian constant 
r1 = sqrt((x+mu).^2 + y.^2 + z.^2);
r2 = sqrt((x + mu -1).^2 + y.^2 + z.^2);
v_squared = vx.^2 + vy.^2 + vz.^2;
C = (0.5*(x.^2 + y.^2) + (1-mu)./r1 + mu./r2 +0.5.*mu.*(1-mu)).*2 - v_squared;
end
function dx = CR3BP_STM_3D(~,xx,mu)

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
function [xf,PHIf,tf, xx, tt]  = propagate(t0,x0,tf,mu,varargin)

    if nargin>4
        evtFlag=varargin{1};
    else
        evtFlag=true;
    end

    tof = tf - t0;

    % Initialize State Transition Matrix at t0
    Phi0 = eye(6);

    % Append to initial conditions the conditions for the STM
    x0Phi0 = [x0; Phi0(:)];
    
    % Perform integration
    options_STM = odeset('reltol', 1e-8, 'abstol', 1e-8,'Events', @(x,y) x_axis_crossing(x,y,evtFlag));
    [tt, xx] = ode78(@(t,x) xyCR3BP_STM_3D(t,x,mu), [0 tof], x0Phi0, options_STM);

    % Extract state vector and State Transition Matrix
    xf = xx(end,1:6)';
    PHIf = reshape(xx(end,7:end),6,6);
    tf = tt(end);

end 
function [value, isterminal, direction] = x_axis_crossing(~,xx,isTerminal)
    value = xx(2);
    isterminal = isTerminal;
    %isterminal=1
    direction = 0;
end
