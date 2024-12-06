%Spacecraft Guidance and Navigation
%Assignment #1
%Author: Mina Baniamein
%% Ex.1 - Guess Generation
clear; clc; close;

%Data
mu=1.21506683*1e-2;
alpha = 0.2*pi;
beta = 1.41;
delta=4;
ti=2;

%Guess generation in the rotating frame
[x0, Ti, Tof]=firstguess_gen(alpha, beta, delta, ti, mu);
Tf = Ti + Tof;
%Change of r.f.: guess in the inertial frame Earth-centered
[X0]=change_frame(x0, Ti, mu);

%Integration
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[tt, xx] = ode113(@(t,x) PBRFBP_rhs(t, x, mu), [Ti Tf], x0, options);

%Plot in the rotating frame
figure(1)
plot(xx(:,1), xx(:,2),'b',...
    'LineWidth',2);
xlabel('x [DU]');
ylabel('y [DU]');
grid on;
real axis;
hold on;
plot(-mu, 0,'k.','MarkerSize',10);
plot(1-mu,0,'k.','MarkerSize',10);
text(-mu, 0,'EARTH','HorizontalAlignment',...
    'center','VerticalAlignment','top');
text(1-mu,0,'MOON','HorizontalAlignment',...
    'center','VerticalAlignment','top');

%Plot in the inertial r.f. earth-centered
disc=length(xx);
XX=zeros(disc,4);
for i=1:length(XX)
    XX(i,:)=change_frame(xx(i,:), tt(i), mu);
end
figure(2)
plot(XX(:,1), XX(:,2),'b', 'LineWidth',2);
xlabel('x [DU]');
ylabel('y [DU]');
grid on;
real axis;
hold on;
plot(0, 0,'k.','MarkerSize',10);
text(0, 0,'EARTH','HorizontalAlignment',...
    'center','VerticalAlignment','top');


%% Ex.2A - Single Shooting without Derivatives
clc; close all;

%Solving NLP
Xopt0=[x0(1:2); x0(3:4); Ti; Tf];
options=optimset('Algorithm','active-set','LargeScale','on',...
    'Display','iter','TolCon',1e-10, 'MaxIter',100,'MaxFunEvals',2000);
[Xopt, ~]=fmincon(@obj_ss, Xopt0, [],[],[],[],[],[], @nonlcons_ss,...
    options);

%Integration
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, xx1] = ode113(@(t,x) PBRFBP_rhs(t, x, mu), [Xopt(5) Xopt(6)],...
    [Xopt(1); Xopt(2); Xopt(3); Xopt(4)], options);

%Plot in the rotating frame
figure(1);
plot(xx1(:,1), xx1(:,2),'b', 'LineWidth',2);
xlabel('x [DU]');
ylabel('y [DU]');
grid on;
real axis;
hold on;
plot(-mu, 0, 'k.','MarkerSize',10);
plot(1-mu,0,'k.','MarkerSize',10);
text(-mu, 0, 'EARTH','HorizontalAlignment',...
    'center','VerticalAlignment','top');
text(1-mu,0,'MOON','HorizontalAlignment',...
    'center','VerticalAlignment','top');

%% Ex.2B - Single Shooting with Derivatives
clc; close all;

%Solving NLP
Xopt0=[x0(1:2); x0(3:4); Ti; Tf];
options=optimset('Algorithm','active-set','LargeScale','on',...
    'Display','iter','TolCon',1e-12,'MaxIter',100,'MaxFunEvals',2000,...
    'GradObj','on','GradConstr','on');
[Xopt, ~]=fmincon(@obj_ss2, Xopt0, [],[],[],[],[],[],...
    @nonlcons_ss2, options);

%Integration
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, xx2] = ode113(@(t,x) PBRFBP_rhs(t, x, mu), [Xopt(5) Xopt(6)],...
    [Xopt(1); Xopt(2); Xopt(3); Xopt(4)], options);

%Plot in the rotating frame of both Single Shooting solutions
figure(1);
plot(xx2(:,1), xx2(:,2),'r', 'LineWidth',2);
xlabel('x [DU]');
ylabel('y [DU]');
ylim([-2.5,3]);
grid on;
hold on;
real axis;
plot(-mu, 0,'k.','MarkerSize',10);
plot(1-mu,0,'k.','MarkerSize',10);
text(-mu, 0, 'EARTH','HorizontalAlignment',...
    'center','VerticalAlignment','top');
text(1-mu,0,'MOON','HorizontalAlignment',...
    'center','VerticalAlignment','top');

%% Functions

function [x0,Ti,ToF]=firstguess_gen(alpha, beta, delta, ti, mu)
%Data
Re=6378;            %Earth radius
hi=167;             %Earth Parking orbit
DU=3.84405000*1e5;  %Distance unit
r0=(Re+hi)/DU;

%Initialization
x0=zeros(4,1); 
v0=beta*sqrt((1-mu)/r0);             %Velocity magnitude (already with DU)

%Guess generation
x0(1)=r0*cos(alpha)-mu;              %State vector
x0(2)=r0*sin(alpha);
x0(3)=-(v0-r0)*sin(alpha);
x0(4)=(v0-r0)*cos(alpha);
Ti=ti;                               %Initial/Final time
ToF=delta;
end

function [X0]=change_frame(x0, Ti, mu)
%Initialization
X0=zeros(4,1);

%Unpack
x=x0(1);
y=x0(2);
xdot=x0(3);
ydot=x0(4);

%Change of r.f.
X0(1)=(x+mu)*cos(Ti)-y*sin(Ti);
X0(2)=(x+mu)*sin(Ti)+y*cos(Ti);
X0(3)=(xdot-y)*cos(Ti)-(ydot+x+mu)*sin(Ti);
X0(4)=(xdot-y)*sin(Ti)+(ydot+x+mu)*cos(Ti);
end

function [dxdt] = PBRFBP_rhs(t, x_vec, mu)
% Initialize
dxdt=zeros(4,1);  %I always initialize the output

x=x_vec(1);
y=x_vec(2);
xdot=x_vec(3);
ydot=x_vec(4);

r1_norm=sqrt((x+mu)^2+y^2);
r2_norm=sqrt((x+mu-1)^2+y^2);

%Additional data for PBRFBP
ms=3.28900541*1e5;
rho=3.88811143*1e2;
ws=-9.25195985*1e-1;

%Composition of a PCRTBP
dxdt(1:2)=[xdot;ydot];
dxdt(3:4)=[2*ydot+x-(1-mu)*(x+mu)/r1_norm^3-mu*(x+mu-1)/r2_norm^3;...
    -2*xdot+y-(1-mu)*y/r1_norm^3-mu*y/r2_norm^3];

%BRFBP through PCRTBP perturbation
r3=sqrt((x-rho*cos(ws*t))^2 + (y-rho*sin(ws*t))^2);
dxdt(3)=dxdt(3)-ms*(x-rho*cos(ws*t))/r3^3-ms*cos(ws*t)/rho^2;
dxdt(4)=dxdt(4)-ms*(y-rho*sin(ws*t))/r3^3-ms*sin(ws*t)/rho^2;

end

function [dxMdt]=STMexact_PBRFBP(t,xM,mu)
%Initialize
dxMdt=zeros(20,1);

x=xM(1);
y=xM(2);
xdot=xM(3);
ydot=xM(4);

M=(reshape(xM(5:20),4,4))';   %From equations to STM

%Additional data for PBRFBP
ms=3.28900541*1e5;
rho=3.88811143*1e2;
ws=-9.25195985*1e-1;

%PCRTBP dynamics
r1_norm=sqrt((x+mu)^2+y^2);
r2_norm=sqrt((x+mu-1)^2+y^2);

dxMdt(1:2)=[xdot;ydot];
dxMdt(3:4)=[2*ydot+x-(1-mu)*(x+mu)/r1_norm^3-mu*(x+mu-1)/r2_norm^3;...
    -2*xdot+y-(1-mu)*y/r1_norm^3-mu*y/r2_norm^3];

%PBRFBP through PCRTBP perturbation
r3=sqrt((x-rho*cos(ws*t))^2 + (y-rho*sin(ws*t))^2);
dxMdt(3)=dxMdt(3)-ms*(x-rho*cos(ws*t))/r3^3-ms*cos(ws*t)/rho^2;
dxMdt(4)=dxMdt(4)-ms*(y-rho*sin(ws*t))/r3^3-ms*sin(ws*t)/rho^2;

%Variational equations of PCRTBP
df1dx=1-(1-mu)/r1_norm^3+3*(1-mu)*(x+mu)^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*(x+mu-1)^2/r2_norm^5;
df1dy=3*(1-mu)*(x+mu)*y/r1_norm^5+3*mu*(x+mu-1)*y/r2_norm^5;
df2dy=1-(1-mu)/r1_norm^3+3*(1-mu)*y^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*y^2/r2_norm^5;

%Variational equations of BRFBP
df1dx=df1dx+ms*(2*(x-rho*cos(ws*t))^2-(y-rho*sin(ws*t))^2)/r3^5;
df1dy=df1dy+3*ms*(x-rho*cos(ws*t))*(y-rho*sin(ws*t))/r3^5;
df2dy=df2dy+ms*(2*(y-rho*sin(ws*t))^2-(x-rho*cos(ws*t))^2)/r3^5;

A=[0, 0, 1, 0;...
    0, 0, 0, 1;...
    df1dx, df1dy, 0, 2;...
    df1dy, df2dy, -2, 0];

dMdt=A*M;
dxMdt(5:8)=dMdt(1,1:4)';
dxMdt(9:12)=dMdt(2,1:4)';
dxMdt(13:16)=dMdt(3,1:4)';
dxMdt(17:20)=dMdt(4,1:4)';
end

function [f] = obj_ss(Xopt)
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
[~, Xf] = ode113(@(t,x) PBRFBP_rhs(t, x, mu), [Ti Tf], [xi; yi; xidot;...
    yidot], options);
Xf=Xf(end,:);

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
end

function [c, ceq] = nonlcons_ss(Xopt)
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
[~, Xf] = ode113(@(t,x) PBRFBP_rhs(t, x, mu), [Ti Tf], [xi; yi; xidot;...
    yidot], options);
Xf=Xf(end, :);

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
end

function [f, df] = obj_ss2(Xopt)
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
[~, xxM] = ode113(@(t,xM) STMexact_PBRFBP(t,xM, mu), [Ti Tf],...
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
    rhs1=PBRFBP_rhs(Ti,[xi; yi; xidot; yidot], mu);
    rhs2=PBRFBP_rhs(Tf,Xf, mu);

    %Components df
    dfdX=[vyi/vi; -vxi/vi; vxi/vi; vyi/vi]+STMf'*...
        [vyf/vf; -vxf/vf; vxf/vf; vyf/vf];
    dfdti=-(STMf*rhs1)'*[vyf/vf; -vxf/vf; vxf/vf; vyf/vf];
    dfdtf=rhs2'*[vyf/vf; -vxf/vf; vxf/vf; vyf/vf];

    df=[dfdX; dfdti; dfdtf];
end
end

function [c, ceq, gradc, gradceq] = nonlcons_ss2(Xopt)
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
[~, xxM] = ode113(@(t,xM) STMexact_PBRFBP(t,xM, mu), [Ti Tf],...
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
    rhs1=PBRFBP_rhs(Ti,[xi; yi; xidot; yidot], mu);
    rhs2=PBRFBP_rhs(Tf,Xf, mu);

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

