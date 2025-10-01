%Spacecraft Guidance and Navigation
%Assignment #2
%Author: Mina Baniamein

%% Point 1
clear ; close all; clc;
%% Kernel unload
cspice_kclear(); %unload kernel section
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));
%% Kernel load
cspice_furnsh('kernels\assignment2.txt'); 
%% LinCov
%Initial data
ri = [-0.011965533749906, -0.017025663128129];
vi = [10.718855256727338, 0.116502348513671];
ti = 1.282800225339865;
tf = 9.595124551366348;  
     

x0 = [ri';vi']; %column vector 4x1 

P0 = [+1.041e-15, +6.026e-17, +5.647e-16, +4.577e-15;  
      +6.026e-17, +4.287e-18, +4.312e-17, +1.855e-16; 
      +5.647e-16, +4.312e-17, +4.432e-16, +1.455e-15;
      +4.577e-15, +1.855e-16, +1.455e-15, +2.822e-14];

% Computation of the PBRFBP mu constant
gm_moon = cspice_bodvrd('MOON', 'GM', 1);
gm_e = cspice_bodvrd('EARTH', 'GM', 1);
mu = gm_moon/(gm_moon + gm_e); 

times = linspace(ti,tf,5);
n = length(times) - 1;

x0_LV=x0;
P0_LV=P0;

mean_statesLV = zeros(4,n);
covariancesLV = zeros(4,4,n);
for i=1:n
    [xf_LV,Pf_LV] = LinCov(x0_LV,times(1),times(i+1),P0_LV,mu);
    mean_statesLV(:,i)   =  xf_LV;
    covariancesLV(:,:,i) =  Pf_LV;
end

mean_statesLV;
covariancesLV;

options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt, xx] = ode113(@(t,x) PBRFBP(t,x,mu), [ti tf], x0_LV, options);

figure
hold on
plot(xx(:,1),xx(:,2),'--')
plot(x0_LV(1),x0_LV(2),'*')
plot(mean_statesLV(1,1),mean_statesLV(2,1),'*')
plot(mean_statesLV(1,2),mean_statesLV(2,2),'*')
plot(mean_statesLV(1,3),mean_statesLV(2,3),'*')
plot(mean_statesLV(1,4),mean_statesLV(2,4),'*')
legend('','t_0','t1','t2','t3','t_f')

figure
hold on
grid on
plotcov(covariancesLV(1:2,1:2,4), mean_statesLV(1:2,4)')
plot(mean_statesLV(1,4),mean_statesLV(2,4),'rx')

%% UT
alpha = 1;
beta = 2;
lambda = n*alpha^2 - n;
times = linspace(ti,tf,5);
n = length(times) - 1;
x0_UT=x0;
P0_UT=P0;

mean_statesUT = zeros(4,n);
covariancesUT = zeros(4,4,n);

for i=1:n
    [y_m,Pf_UT] = UT(x0_UT,times(1),times(i+1),P0_UT,alpha,beta,mu);
    mean_statesUT(:,i)   =  y_m;
    covariancesUT(:,:,i) =  Pf_UT;
end

mean_statesUT;
covariancesUT;

figure
hold on
plot(xx(:,1),xx(:,2),'--')
plot(x0_LV(1),x0_LV(2),'*')
plot(mean_statesUT(1,1),mean_statesUT(2,1),'*')
plot(mean_statesUT(1,2),mean_statesUT(2,2),'*')
plot(mean_statesUT(1,3),mean_statesUT(2,3),'*')
plot(mean_statesUT(1,4),mean_statesUT(2,4),'*')
legend('','t_0','t1','t2','t3','t_f')

figure
hold on
grid on
error_ellipse(covariancesLV(1:2,1:2,4), mean_statesLV(1:2,4), 'b')
error_ellipse(covariancesUT(1:2,1:2,4), mean_statesUT(1:2,4), 'r')
plot(mean_statesLV(1,4),mean_statesLV(2,4),'g*')
plot(mean_statesUT(1,4),mean_statesUT(2,4),'mx')
legend('LV ellipse', 'UT ellipse', 'LV mean', 'UT mean')


%% Point 2
% Given data
x0; P0;                  % Covariance matrix
N = 1000;                % Number of Monte Carlo samples

rng('default'); 
samples = mvnrnd(x0, P0, N);  %1000x4 with initial states as rows
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
propagated_samples = zeros(N,4,4);

n = length(x0);
mean_statesMC = zeros(4,n);
covariancesMC = zeros(4,4,n);
ti = times(1);

for s=1:n
    tf = times(s+1);
    for i=1:N
        [~, xx] = ode113(@(t,x) PBRFBP(t,x,mu), [ti tf], samples(i,:)', options); %required column initial states
        propagated_samples(i,:,s) = xx(end,:);
    end
    mean_statesMC(:,s) = mean(propagated_samples(:,:,s));
    covariancesMC(:,:,s) = cov(propagated_samples(:,:,s));
end
%% MC propagation and ellipses

figure 
hold on
grid on
plot(xx(:,1),xx(:,2),'--')
plot(x0_LV(1),x0_LV(2),'*')
plot(propagated_samples(:, 1, 1), propagated_samples(:, 2, 1), '*');
plot(propagated_samples(:, 1, 2), propagated_samples(:, 2, 2), '*');
plot(propagated_samples(:, 1, 3), propagated_samples(:, 2, 3), '*');
plot(propagated_samples(:, 1, 4), propagated_samples(:, 2, 4), '*');
legend('','t_0','t1','t2','t3','t_f')
title('Montecarlo propagation')


mean_statesMC;
covariancesMC;

figure
hold on
grid on
error_ellipse(covariancesMC(1:2,1:2,4), mean_statesMC(1:2,4), 'g')
error_ellipse(covariancesLV(1:2,1:2,4), mean_statesLV(1:2,4), 'b')
error_ellipse(covariancesUT(1:2,1:2,4), mean_statesUT(1:2,4), 'r')
plot(mean_statesMC(1,4),mean_statesMC(2,4),'rx') 
plot(mean_statesLV(1,4),mean_statesLV(2,4),'g*')
plot(mean_statesUT(1,4),mean_statesUT(2,4),'mx')
legend('MC ellipse','LV ellipse', 'UT ellipse','MC mean', 'LV mean', 'UT mean')


%% Pr eigs plot
figure
hold on
grid on
xlabel('X Position');
ylabel('Y Position');
title('Eigenvalues PR Evolution');

% Plot initial point for all 3 methods
plot(times(1), 3*sqrt(max(eig(P0(1:2,1:2)))), 'bo', 'DisplayName', 'yLV Initial');
plot(times(1), 3*sqrt(max(eig(P0(1:2,1:2)))), 'gx', 'DisplayName', 'yUT Initial');
plot(times(1), 3*sqrt(max(eig(P0(1:2,1:2)))), 'rs', 'DisplayName', 'yMC Initial');


for i = 1:n
    yLV = 3*sqrt(max(eig(covariancesLV(1:2,1:2,i))));
    yUT = 3*sqrt(max(eig(covariancesUT(1:2,1:2,i))));
    yMC = 3*sqrt(max(eig(covariancesMC(1:2,1:2,i))));

    % Plot each dataset
    plot(times(i+1), yUT, 'gx', 'DisplayName', 'yUT');
    plot(times(i+1), yMC, 'rs', 'DisplayName', 'yMC');
    plot(times(i+1), yLV, 'bo', 'DisplayName', 'yLV');
end

legend('yLV', 'yUT', 'yMC');

%% Pv eigs plot
figure
hold on
grid on
xlabel('X Position');
ylabel('Y Position');
title('Eigenvalues PV Evolution');

% Plot initial point for all 3 methods
plot(times(1), 3*sqrt(max(eig(P0(3:4,3:4)))), 'bo', 'DisplayName', 'yLV Initial');
plot(times(1), 3*sqrt(max(eig(P0(3:4,3:4)))), 'gx', 'DisplayName', 'yUT Initial');
plot(times(1), 3*sqrt(max(eig(P0(3:4,3:4)))), 'rs', 'DisplayName', 'yMC Initial');


for i = 1:n
    yLV = 3*sqrt(max(eig(covariancesLV(3:4,3:4,i))));
    yUT = 3*sqrt(max(eig(covariancesUT(3:4,3:4,i))));
    yMC = 3*sqrt(max(eig(covariancesMC(3:4,3:4,i))));

    % Plot each dataset
    plot(times(i+1), yUT, 'gx', 'DisplayName', 'yUT');
    plot(times(i+1), yMC, 'rs', 'DisplayName', 'yMC');
    plot(times(i+1), yLV, 'bo', 'DisplayName', 'yLV');
end

legend('yLV', 'yUT', 'yMC');


%% QQ plots
figure
hold on
grid on
qqplot(propagated_samples(:,1,end))
title('QQplot x');

figure
hold on
grid on
qqplot(propagated_samples(:,2,end))
title('QQplot y');

figure
hold on
grid on
qqplot(propagated_samples(:,3,end))
title('QQplot v_x');

figure
hold on
grid on
qqplot(propagated_samples(:,4,end))
title('QQplot v_y');

%% Functions
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

function [xf,Pf] = LinCov(x0,ti,tf,P0,mu)
% DESCRIPTION
%   This function propagates the state and covariance matrix using linearized 
%   equations of motion based on the Planar Bicircular Restricted Four-Body Problem (PBRFBP)
%   and computes the final state and covariance at a specified final
%   time, based on the LinCov method
%
% INPUT
%   x0  [4x1]   Initial adimensionalized state vector (position and velocity)              [-]
%   ti  [1x1]   Initial time                                                               [s]
%   tf  [1x1]   Final time                                                                 [s]
%   P0  [4x4]   Initial state covariance matrix                                            [-]
%   mu  [1x1]   Adimensionalized mass parameter                                            [-]
%
% OUTPUT
%   xf  [4x1]   Final propagated mean state vector                                         [-]
%   Pf  [4x4]   Final propagated state covariance matrix                                   [-]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

phi0 = reshape(eye(4),16,1);
xx0  = [x0;phi0]; 
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[~, xx] = ode113(@(t,x) PBRFBP_STM(t,x,mu), [ti tf], xx0, options);
xf = xx(end,1:4)'; %mean final propagated state
PHI = reshape(xx(end,5:end), 4,4);
Pf = PHI*P0*PHI'; %final covariance
end

function [y_m,Py] = UT(x0,ti,tf,P0,alpha,beta,mu)
% DESCRIPTION
%   This function performs uncertainty propagation using the Unscented Transform (UT)
%   for the Planar Bicircular Restricted Four-Body Problem (PBRFBP). It generates sigma 
%   points, propagates them through the nonlinear dynamics, and computes the propagated 
%   mean state and covariance matrix.
%
% INPUT
%   x0     [4x1]   Initial adimensionalized state vector (position and velocity)             [-]
%   ti     [1x1]   Initial time                                                              [s]
%   tf     [1x1]   Final time                                                                [s]
%   P0     [4x4]   Initial state covariance matrix                                           [-]
%   alpha  [1x1]   Tuning parameter                                                          [-]
%   beta   [1x1]   Tuning parameter                                                          [-]
%   mu     [1x1]   Adimensionalized mass parameter                                           [-]
%
% OUTPUT
%   y_m    [4x1]   Propagated mean state vector                                              [-]
%   Py     [4x4]   Propagated state covariance matrix                                        [-]
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

n = length(x0);
lambda = n*alpha^2 - n;
P_scaled = sqrtm((n+lambda)*P0);

%chi points computation
chi = zeros(4,2*n +1); %inizialization of a 4x9 matrix to store all chi column vectors
chi(:,1) = x0; % 4x1
for i=1:n
    chi(:,1+i) = x0 + P_scaled(:,i);
    chi(:,1+n+i) = x0 - P_scaled(:,i);
end

%chi points propagation
Y = zeros(4,2*n +1); % 4x9 matrix to store propagations
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
for i=1:(2*n+1)
    [~, xx] = ode113(@(t,x) PBRFBP(t,x,mu), [ti tf], chi(:,i), options); 
    Y(:,i) = [xx(end,:)]';
end

% weights
Wm0 = lambda/(n+lambda);
Wc0 = lambda/(n+lambda) + 1-alpha^2+beta;
Wi  = 1/(2*(n+lambda));

Wm = [Wm0, Wi.*ones(1,2*n)];
Wc = [Wc0, Wi.*ones(1,2*n)];

y_m = zeros(4,1); %inizialization of mean state
Py = zeros(4,4);  %inizialization of covariance
for i=1:(2*n + 1)
    y_m = y_m + Wm(1,i).*Y(:,i);
end
for i=1:(2*n + 1)
    Py  = Py + Wc(1,i).*((Y(:,i)-y_m)*( Y(:,i)-y_m )');
end
end

function error_ellipse(P, mean, color)
% DESCRIPTION
%   This function plots a 2D error ellipse representing confidence region 
%   of a covariance matrix, centered at a given mean state. The ellipse is computed 
%   from the eigenvalues and eigenvectors of the covariance matrix.
%
% INPUT
%   P      [2x2]   Covariance matrix for the two variables being plotted                    [-]
%   mean   [1x2]   Mean state vector (center of the ellipse)                                [-]
%   color  [1x1]   Line color/style specification for the plot (e.g., 'r', 'b--')           [-]
%
% OUTPUT
%   -              This function produces a plot and does not return any output.
%
% AUTHOR
%   Mina Baniamein, A.Y. 2024/25, MATLAB

[V, D] = eig(P);
t = linspace(0, 2*pi, 100);
ellip = [cos(t); sin(t)]' * sqrt(D) * V';
plot(ellip(:, 1) + mean(1), ellip(:, 2) + mean(2), color);
end



