%% data

% s/c data
mu_e = astroConstants(13);

m=14;
w=10*10^-2;
d=20*10^-2;
h=30*10^-2;

I_x=0.050566;
I_y=0.12641;
I_z=0.5517;
J = [I_x, 0,0; 0, I_y, 0; 0,0,I_z];

tf = 1000;
A0 = [0 1 0; 0 0 1; 1, 0, 0];

w0_y = 0.7;
w0_z = 0.2;
w0_x = 0.2;
w0 = [w0_x; w0_y; w0_z];

% Verifying Inertia moments
h_v = J * w0;
h = norm(h_v);
T = w0' * h_v;

I_vec = [I_x, I_y, I_z];

I_max = max(I_vec);
I_min = min(I_vec);
for i = 1:3
    if I_min < I_vec(i) && I_max > I_vec(i)
        I_int = I_vec(i);
    end
end

if I_max > ((h^2) / (2*T)) && I_min < ((h^2) / (2*T))
    disp('Inertia moment values are acceptable')
else
    disp('WARNING: Inertia moment values not acceptable')
end

% Orbital data
Radius_E = astroConstants(23);
a_sc = (((100 * 60)/(2*pi))^2 * mu_e)^(1/3);
e_sc = 0;
i_sc = 90 * pi / 180;
n_sc = sqrt(mu_e / (a_sc^3));

R_S2E = astroConstants(2);

epsilon = 23.45 * pi / 180;

w_LN = [0;0;n_sc];

n_sun = 2*pi / (astroConstants(32) * 24 * 60 * 60);


% Solar radiation pressure data

Ndir_b1 = [1;0;0];
Ndir_b2 = [0;1;0];
Ndir_b3 = [-1;0;0];
Ndir_b4 = [0;-1;0];
Ndir_b5 = [0;0;1];
Ndir_b6 = [0;0;-1];

Ndir_p1 = [1;0;0];
Ndir_p2 = [-1;0;0];
Ndir_p3 = [1;0;0];
Ndir_p4 = [-1;0;0];

rhoS_b = 0.5;
rhoS_p = 0.1;

rhoD = 0.1;

A_b1 = 6e-2;
A_b2 = 6e-2;
A_b3 = 6e-2;
A_b4 = 6e-2;

A_b5 = 4e-2;
A_b6 = 4e-2;

A_p = 12e-2;

F_e = 1358; % w/m^2
c = astroConstants(5) * 10^3; % km/s

P = F_e / c;

w = 0.2;
d = 0.2;
h = 0.3;
l_panel = 0.6;


rf_b1 = [w/2;0;0];
rf_b2 = [0;d/2;0];
rf_b3 = [-w/2;0;0];
rf_b4 = [0;-d/2;0];
rf_b5 = [0;0;h/2];
rf_b6 = [0;0;-h/2];

rf_p1 = [0;0;h/2 + l_panel/2];
rf_p2 = [0;0;h/2 + l_panel/2];
rf_p3 = [0;0;-(h/2 + l_panel/2)];
rf_p4 = [0;0;-(h/2 + l_panel/2)];


% Magnetic torque data

j_b = [0.0056; 0.0056; 0.0056];

Radius_E = astroConstants(23);

g_10 = -29404.8;
g_11 = -1450.9;
h_11 = 4652.5;

H0 = sqrt((g_10^2) + (g_11^2) + (h_11^2));

w_earth = (15.04 * pi / 180) / 3600;
%w_earth = 2 * pi / (24*60*60);

T_orbit = 2 * pi / n_sc;

%% 2Bp data

a_2bp = a_sc;
e_2bp = e_sc;
i_2bp = i_sc;
om = 0;
OM = pi/2;
theta0 = 0;

[r_0,v_0] = kep2car (a_2bp,e_2bp,i_2bp,OM,om,theta0,mu_e);

mu_E = astroConstants(13);

 % [km/s]
y0 = [ r_0; v_0 ];
% Set time span
a = 1/( 2/norm(r_0) - dot(v_0,v_0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a_sc^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 5*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );

%% control logic

% R = [0 -d/2 h/2 0 -d/2 h/2 0 d/2 -h/2 d/2 0 -d/2; -h/2 0 -w/2 -h/2 0 -w/2 -h/2 0 w/2 0 -h/2 -w/2; w/2 -w/2 0 -w/2 w/2 0 w/2 -w/2 0 -w/2 w/2 0];
Kp = 0.01;
Kd = 0.05;
Tu = 0;

RR = [0.0851 0.0851 -0.0851 -0.0851; -0.0525 0.0525 0.0525 -0.0525; 0.0326 -0.0326 0.0326 -0.0326];
wnull = null(RR, 'rational');

ARW = 4.365 *10^-5;
RRW = 1.455 *10^-6;
A_e_gyro = [0 -(0.04*pi)/180 (0.04*pi)/180; (0.04*pi)/180 0 -(0.04*pi)/180; -(0.04*pi)/180 (0.04*pi)/180 0];
A_e = eye(3)-A_e_gyro;

k1_slewmanoeuvre = 1;
k2_slewmanoeuvre = 0.05;
K_detumbling = 0.5;
