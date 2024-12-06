%Spacecraft Guidance and Navigation
%Assignment #1
%Author: Mina Baniamein

%% Exercise 1
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
gm_moon = cspice_bodvrd('MOON', 'GM', 1);
gm_e = cspice_bodvrd('EARTH', 'GM', 1);
mu = gm_moon/(gm_moon + gm_e);

ri = (Re + hi)/DU;
rf = (R_moon + hf)/DU;
r0 = ri;
v0 = beta*sqrt((1-mu)/r0);


x0 = r0*cos(alpha) - mu;
y0 = r0*sin(alpha);
vx0 = -(v0 -r0)*sin(alpha);
vy0 = (v0 - r0)*cos(alpha);

xx0 = [x0, y0, vx0, vy0]';

options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt, xx] = ode113(@(t,x) PBRFBP(t,x,mu), [2 6], xx0, options);

figure('Name','Rotating Frame') %plottare Luna
plot(xx(:,1),xx(:,2))

n = size(xx,1);
xx_ECI = zeros(n,4);

for i=1:n
    ECI = rot_to_ECI(xx(i,:), mu, tt(i));
    xx_ECI(i,1) = ECI(1);
    xx_ECI(i,2) = ECI(2);
    xx_ECI(i,3) = ECI(3);
    xx_ECI(i,4) = ECI(4);
end

figure('Name','ECI') %plottare Luna
plot(xx_ECI(:,1),xx_ECI(:,2))

%check if frame change is correct
err = zeros(n,1);
for i=1:n
    err(i) = abs(norm([xx(i,1)+mu,xx(i,2)]) - norm([xx_ECI(i,1), xx_ECI(i,2)]));
end



figure('Name','err')
plot(linspace(0,n,n), err)
%% Point 2.a
xx0 = [xx0; 2; 6];

options=optimset('Algorithm','active-set','LargeScale','on','Display','iter','TolCon',1e-10 'MaxIter',500,'MaxFunEvals',2000);  %MODIFICARE OPZIONI PRIMA DELLA CONSEGNA
[xx0_optimal, ~] = fmincon(@deltavfun,xx0,[],[],[],[],[],[],@nonlcon,options) %NLP solving
 
delta_v = deltavfun(xx0_optimal)

options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt, xx_optimal] = ode113(@(t,x) PBRFBP(t,x,mu), [xx0_optimal(5) xx0_optimal(6)], xx0_optimal(1:4), options); %optimal propagation

figure('Name','OPTIMAL - NO STM  - Rotating Frame') %add on plot Moon
plot(xx_optimal(:,1),xx_optimal(:,2))

%% Point 2.b
Xopt0=[xx0(1:2); xx0(3:4); 2; 6];
options=optimset('Algorithm','active-set','LargeScale','on',...
    'Display','iter','TolCon',1e-10,'MaxIter',100,'MaxFunEvals',2000,...
    'GradObj','on','GradConstr','on');
[Xopt, ~]=fmincon(@deltavfun_grad, Xopt0, [],[],[],[],[],[],...
    @nonlcons_grad, options);
options=optimset('Algorithm','active-set','LargeScale','on','Display','iter','TolCon',1e-12, 'MaxIter',100,'MaxFunEvals',10000,'GradObj','on','GradConstr','on');  %MODIFICARE OPZIONI PRIMA DELLA CONSEGNA

%xx0_optimal_grad = fmincon(@deltavfun_grad,xx0,[],[],[],[],[],[],@nonlcons_grad,options)

delta_v = deltavfun(Xopt)


%% Functions

function [dxdt] = PBRFBP(t,xx, mu)

    
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

    A=[0, 0, 1, 0;           %prima di jacobians definire A come zeros e cambiare elementi 
       0, 0, 0, 1;
       df1dx, df1dy, 0, 2;
       df1dy, df2dy, -2, 0];

    dPHI = A*PHI;

    dxdt(5:end) = reshape(dPHI,16,1);
end
function [xx_ECI] = rot_to_ECI(xx_rot, mu, t) %usare un for CASO 2D

xx_ECI = zeros(1,4);

xx_ECI(1) = (xx_rot(1) + mu)*cos(t) - xx_rot(2)*sin(t);
xx_ECI(2) = (xx_rot(1) + mu)*sin(t) + xx_rot(2)*cos(t);
xx_ECI(3) = (xx_rot(3) - xx_rot(2))*cos(t) - (xx_rot(4) + xx_rot(1) + mu)*sin(t);
xx_ECI(4) = (xx_rot(3) - xx_rot(2))*sin(t) + (xx_rot(4) + xx_rot(1) + mu)*cos(t);

end
function delta_v = deltavfun(xx_i)

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
[tt, xx] = ode113(@(t,x) PBRFBP(t,x,mu), [ti tf], [xx_i(1);xx_i(2);xx_i(3);xx_i(4)], options);
xx_f = xx(end, :);

delta_v = (sqrt( (xx_i(3) -xx_i(2))^2 + (xx_i(4) + xx_i(1) + mu)^2) - sqrt((1-mu)/ri)) + abs(sqrt((xx_f(3) - xx_f(2))^2 + (xx_f(4) + xx_f(1) + mu - 1)^2) - sqrt(mu/rf)); 

end
function [c, ceq] = nonlcon(xx_i)

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
[tt, xx] = ode113(@(t,x) PBRFBP(t,x,mu), [ti tf], [xx_i(1);xx_i(2);xx_i(3);xx_i(4)], options);
xx_f = xx(end, :);

ceq(1) = (xx_i(1) + mu)^2 + xx_i(2)^2 - ri^2;
ceq(2) = (xx_i(1)+mu)*(xx_i(3) - xx_i(2)) + xx_i(2)*(xx_i(4) + xx_i(1) + mu);
ceq(3) = (xx_f(1) + mu -1)^2 + xx_f(2)^2 - rf^2;
ceq(4) = (xx_f(1)+mu-1)*(xx_f(3) - xx_f(2)) + xx_f(2)*(xx_f(4) + xx_f(1) + mu-1);

c = [];
end
function [f, df] = deltavfun_grad(Xopt) 
%Data
mu=1.21506683*1e-2;

%Upacking state
xi=Xopt(1);
yi=Xopt(2);
xidot=Xopt(3);
yidot=Xopt(4);
Ti=Xopt(5);
Tf=Xopt(6);

%Earth orbit
Re=6378;            %Earth radius
hi=167;             %Earth Parking orbit
DU=3.84405000*1e5;  %Distance unit
ri=(Re+hi)/DU;

%Integration
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, xxM] = ode113(@(t,xM) PBRFBP_STM(t,xM, mu), [Ti Tf],...
        [[xi; yi; xidot; yidot]; reshape(eye(4),[],1)], options);
Xf=xxM(end,1:4);
STMf=(reshape(xxM(end,5:20),4,4))';   %From equations to STM

%Final state
xf=Xf(1);
yf=Xf(2);
xfdot=Xf(3);
yfdot=Xf(4);

%Moon orbit
Rm=1738;            %Earth radius
hf=100;             %Earth Parking orbit
rf=(Rm+hf)/DU;

%OF
dV1=abs(sqrt((xidot-yi)^2+(yidot+xi+mu)^2)-sqrt((1-mu)/ri));
dV2=abs(sqrt((xfdot-yf)^2+(yfdot+xf+mu-1)^2)-sqrt(mu/rf));
f=dV1+dV2;

%DOF
if nargout>1
    
    %Velocity at Ti
    vi=sqrt((xidot-yi)^2+(yidot+xi+mu)^2);
    vxi=xidot-yi; 
    vyi=yidot+xi+mu;

    %Velocity at Tf
    vf=sqrt((xfdot-yf)^2+(yfdot+xf+mu-1)^2);
    vxf=xfdot-yf; 
    vyf=yfdot+xf+mu-1;
    
    %Accelerations at Ti, Tf in 2D
    rhs1=PBRFBP(Ti,[xi; yi; xidot; yidot], mu);
    rhs2=PBRFBP(Tf,Xf, mu);

    %Components df
    dfdX=[vyi/vi; -vxi/vi; vxi/vi; vyi/vi]+STMf'*...
        [vyf/vf; -vxf/vf; vxf/vf; vyf/vf];
    dfdti=-(STMf*rhs1)'*[vyf/vf; -vxf/vf; vxf/vf; vyf/vf];
    dfdtf=rhs2'*[vyf/vf; -vxf/vf; vxf/vf; vyf/vf];

    df=[dfdX; dfdti; dfdtf];
end
end
function [c, ceq, gradc, gradceq] = nonlcons_grad(Xopt)
%Data
mu=1.21506683*1e-2;

%Upacking state
xi=Xopt(1);
yi=Xopt(2);
xidot=Xopt(3);
yidot=Xopt(4);
Ti=Xopt(5);
Tf=Xopt(6);

%Earth orbit
Re=6378;            %Earth radius
hi=167;             %Earth Parking orbit
DU=3.84405000*1e5;  %Distance unit
ri=(Re+hi)/DU;

%Integration
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, xxM] = ode113(@(t,xM) PBRFBP_STM(t,xM, mu), [Ti Tf],...
        [[xi; yi; xidot; yidot]; reshape(eye(4),[],1)], options);
Xf=xxM(end,1:6);
STMf=(reshape(xxM(end,5:20),4,4))';   %From equations to STM

%Final state
xf=Xf(1);
yf=Xf(2);
xfdot=Xf(3);
yfdot=Xf(4);

%Moon orbit
Rm=1738;            %Earth radius
hf=100;             %Earth Parking orbit
rf=(Rm+hf)/DU;

%C
c=[];

%Ceq
ceq=[(xi+mu)^2+yi^2-ri^2;...
    (xi+mu)*(xidot-yi)+yi*(yidot+xi+mu);...
    (xf+mu-1)^2+yf^2-rf^2;...
    (xf+mu-1)*(xfdot-yf)+yf*(yfdot+xf+mu-1)];

if nargout>2
    %Gradient inequality constraints
    gradc=[];

    %Accelerations at Ti, Tf in 2D
    rhs1=PBRFBP(Ti,[xi; yi; xidot; yidot], mu);
    rhs2=PBRFBP(Tf,Xf, mu);

    %Gradient equality constraints
    gradceq_12=[2*(xi+mu), 2*yi, 0, 0, 0, 0;...
        xidot, yidot, xi+mu, yi, 0, 0];

    gradceq_3_state=2*[xf+mu-1, yf, 0, 0]*STMf;
    gradceq_3_time=2*[xf+mu-1, yf, 0, 0]*[-STMf*rhs1, rhs2];
    gradceq_3=[gradceq_3_state, gradceq_3_time];

    gradceq_4_state=[xfdot, yfdot, xf+mu-1, yf]*STMf;
    gradceq_4_time=[xfdot, yfdot, xf+mu-1, yf]*[-STMf*rhs1, rhs2];
    gradceq_4=[gradceq_4_state, gradceq_4_time];

    gradceq=[gradceq_12; gradceq_3; gradceq_4];
    gradceq=gradceq';
end
end


















