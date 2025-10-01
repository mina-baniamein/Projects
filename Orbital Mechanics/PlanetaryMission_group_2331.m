%% Assignment 2 group 2331

m=3;
k=14;
K=60; % number of orbits
AM=0.0254;
Cd=2.1;
J2=astroConstants(9);
i=deg2rad(54.9556); % [rad]
e=0.5456;
a=15051; % [km]
w=deg2rad(139.9687); % ---------------> given by us
RAAN=deg2rad(263.5589); % ---------------> given by us
nu=deg2rad(0); % ---------------> given by us

t0=0;
theta0=0; 
we=(deg2rad(15.04))/3600;
wevector=[0;0;we];
mu=astroConstants(13);
Re=astroConstants(23);

T = 2*pi*sqrt(a^3/mu); % Orbital period [1/s]
numerore=T/3600; % number of hours
tspan=linspace(0,K*T,100000);

[r0,v0] = kep2cartmio(a,e,i,w,nu,RAAN,mu);
y0=[r0;v0];

% REPEATING GROUND TRACK FOR AN UNPERTURBATED ORBIT
[T1,a1] = repeatinggroundtrack(mu,we,k,m); % here we compute the new a1 and T1
[r1,v1] = kep2cartmio(a1,e,i,w,nu,RAAN,mu); % 
y1=[r1;v1];

fprintf('Semimajoraxis for repeating ground track: %.3f km \n',a1)
tspan2=linspace(0,K*T1,100000); % the new vector of times for the new time period of the orbit

% UNPERTURBATED ORBIT PROPAGATION
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[t, y]=ode113(@(t,y) odefun(t,y,mu),tspan,y0,options); % risolutore equazione differenziale, in cui abbiamo uuna funzione
r_t=y(:,1:3);
r_tf=r_t(end,:);

% GROUND TRACK
[alpha, delta, lon, lat] = groundTrack(t0,tspan,y0,theta0,we,mu);
[alpha1, delta1, lon1, lat1] = groundTrack(t0,tspan2,y1,theta0,we,mu);

% UNPERTURBATED ORBIT PLOT
figure()
plotearth() % PLOTTING THE EARTH
hold on
plot3(y(:,1),y(:,2),y(:,3),'b-','Linewidth',2)
scatter3(r0(1),r0(2),r0(3),100,'green','square','filled')
scatter3(r_tf(1),r_tf(2),r_tf(3),'blue','square','filled')
xlabel('x [Km]',Interpreter='latex');ylabel('y [km]',Interpreter='latex');zlabel('z [km]',Interpreter='latex');
title('Unperturbated orbit',[num2str(K*T/3600) ' hours | ' num2str(K) ' Periods'],Interpreter='latex')
legend('Earth','S/C','Initial point','Final point',Interpreter='latex')
axis equal;
grid on

% GROUND TRACK WITHOUT REPETITIONS AND WITH REPETITIONS
figure()
plot(lon,lat,'LineStyle','none','Marker','.');
hold on
scatter(lon(1),lat(1),500,'green','square','filled')
scatter(lon(end),lat(end),500,'blue','square','filled')
xlim([-180 180]);ylim([-90 90]);
xlabel('Longitude [deg]',Interpreter='latex');ylabel('Latitude [deg]',Interpreter='latex')
plot(lon1,lat1,'LineStyle','none','Marker','.');
scatter(lon1(end),lat1(end),500,'blue','square','filled')
terra=imread('earthh.jpg');
legend('Original','initial position','end position','Repeating','end position (repeating)',Interpreter='latex');
title('Ground track and repeating ground track for an unperturbated orbit',Interpreter='latex')
hold on
EEEARTH=image(xlim,ylim,flipdim(terra,1));
uistack(EEEARTH,'bottom')
grid on
hold off

% PERTURBATED ORBIT PROPAGATION
[tt, yy]=ode113(@(t,y) odefunperturb(t,y,mu,Re,J2,AM,Cd,we,'cart'),tspan,y0,options);
r_t1=yy(:,1:3);
r_tf1=r_t1(end,:);

% PERURBATED GROUND TRACK
[alpha2, delta2, lon2, lat2] = perturbatedgroundTrack(t0,tspan,y0,theta0,we,mu,AM,J2,Cd,Re);
[alpha3, delta3, lon3, lat3] = perturbatedgroundTrack(t0,tspan2,y1,theta0,we,mu,AM,J2,Cd,Re);

nn=100000;
% PERTURBATED ORBIT PLOT
figure()
plotearth()
hold on
%plot3(yy(:,1),yy(:,2),yy(:,3),'b-') %orbita perturbata
patch([yy(:,1) nan(nn,1)],[yy(:,2) nan(nn,1)],[yy(:,3) nan(nn,1)],[t/3600 nan(nn,1)],'FaceColor','none','EdgeColor','interp')
c = colorbar;
c.Label.String = 'time [hours]';
scatter3(r0(1),r0(2),r0(3),100,'green','square','filled')
scatter3(r_tf1(1),r_tf1(2),r_tf1(3),'blue','square','filled')
xlabel('x [Km]',Interpreter='latex');ylabel('y [km]',Interpreter='latex');zlabel('z [km]',Interpreter='latex');
legend('Earth','Perturbated orbit',Interpreter='latex')
title('Perturbated orbit',[num2str(K*T/3600) ' hours | ' num2str(K) ' Periods'],Interpreter='latex')
axis equal;
grid on
hold off

% GROUND TRACK WITHOUT REPETITIONS AND WITH REPETITIONS FOR A PERTURBATED
% ORBIT

figure()
plot(lon2,lat2,'LineStyle','none','Marker','.');
hold on
scatter(lon2(1),lat2(1),500,'green','square','filled')
scatter(lon2(end),lat2(end),500,'blue','square','filled')
xlim([-180 180]);ylim([-90 90]);
xlabel('Longitude [deg]',Interpreter='latex');ylabel('Latitude [deg]',Interpreter='latex')
plot(lon3,lat3,'LineStyle','none','Marker','.');
scatter(lon3(end),lat3(end),500,'blue','square','filled')
terra=imread('earthh.jpg');
legend('Original','initial position','end position','Repeating','end position (repeating)',Interpreter='latex');
title('Ground track and repeating ground track for a perturbated orbit',Interpreter='latex')
hold on
EEEARTH=image(xlim,ylim,flipdim(terra,1));
uistack(EEEARTH,'bottom')
grid on
hold off

%%
% Modeling perturbations with Gauss equations
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
%angular velocity of Earth around its axis:
we = deg2rad(15.04)/3600; %[rad/s]
%planetary constant of the Earth:
mu_e = astroConstants(13); %[km^3/s^2]
%Radius of the Earth:
Re = astroConstants(23); %[km]
a0 = 15051; %semi-major axis [km]
e0 = 0.5456;  %eccentricity
i0 = deg2rad(54.9556); %inclination [rad]
% I have to choose the W,w,f0 parameteres:
OMEGA0 = deg2rad(263.5589); %[rad]
omega0 = deg2rad(139.9687); %[rad]
n=10000;
f0 = deg2rad(0);   %[rad]
T = 2*pi*sqrt( a0^3/mu_e); K = 60;
tspan = linspace(0,K*T,n);
   
% SRP + J2
J2 = astroConstants(9);
AU = astroConstants(2);
Cd = 2.1; % drag coefficient
AM = 0.0254; % spacecraft area-to-mass ratio relative to the air [m^2/kg]

% Propagation of keplerian elements

y0_kepGauss = [a0,e0,i0,OMEGA0,omega0,f0];
[r0, v0] = kep2cartmio(a0,e0,i0,omega0,f0,OMEGA0,mu_e); y0_cart = [r0;v0];
tic
[~,y_cart] = ode113(@(t,y) odefunperturb(t,y, mu_e, Re, J2,AM,Cd,we,'cart'), tspan, y0_cart, options);
tempocomputCart=toc;
tic
[t,y_kepGauss] = ode113(@(t,y) odefunperturb(t,y, mu_e, Re, J2,AM,Cd,we,'RSW') , tspan, y0_kepGauss, options);
tempocomputGauss=toc;
y_kep = [];
for j = 1:length(t)
    [a,e,i,Omega,omega,f]= cart2kep(y_cart(j,1:3),y_cart(j,4:6),mu_e,'rad');
    y_kep(j,:) = [a,e,i,Omega,omega,f];
end
A = movmean(y_kep(:,1), 400);
E = movmean(y_kep(:,2), 400);
I = movmean(rad2deg(y_kep(:,3)), 400);
O = movmean(rad2deg(y_kep(:,4)), 400);
o = movmean(rad2deg(y_kep(:,5)), 400);
F = movmean(rad2deg(unwrap(y_kep(:,6))), 400);

fprintf('\n');
fprintf('Computational time with Cartesian propagation: %.7f sec\n',tempocomputCart);
fprintf('Computational time with Gauss propagation: %.7f sec\n',tempocomputGauss);
fprintf('\n');

figure(7)  %Semi-major axis 
S_major_axis = tiledlayout(2,2);  
title(S_major_axis,'Semi-major axis',Interpreter='latex') 
 
nexttile 
plot(t/T,y_kepGauss(:,1),'b-')
title([num2str(K) ' Periods'],Interpreter='latex')
xlabel('t [T]',Interpreter='latex');ylabel('a [km]',Interpreter='latex')
hold on 
plot(t/T,y_kep(:,1),'r-') 
xlabel('t [T]',Interpreter='latex');ylabel('a [km]',Interpreter='latex')
plot(t/T, A,'-','Color','[1, 0.84, 0]',LineWidth=2)
xlabel('t [T]',Interpreter='latex');ylabel('a [km]',Interpreter='latex')
legend('Gauss equations','Cartesian','Secular (filtered)') 
xlim([0 K]) 
grid on; 
hold off 
 
nexttile  
plot(t(1:n/10)/T,y_kepGauss(1:(n/10),1),'b-')
title([num2str(K/10) ' Periods'],Interpreter='latex')
xlabel('t [T]',Interpreter='latex');ylabel('a [km]',Interpreter='latex')
hold on 
plot(t(1:n/10)/T,y_kep(1:(n/10),1),'r-')
xlabel('t [T]',Interpreter='latex');ylabel('a [km]',Interpreter='latex')
plot(t(1:n/10)/T, A(1:(n/10)),'-','Color','[1, 0.84, 0]',LineWidth=2)
xlabel('t [T]',Interpreter='latex');ylabel('a [km]',Interpreter='latex')
legend('Gauss equations','Cartesian','Secular (filtered)') 
xlim([0 K/10]) 
grid on; 
hold off 
 
nexttile(3,[1 2]) 
semilogy(t/T,abs(y_kep(:,1)-y_kepGauss(:,1))/a0,'b-') 
title('Relative error',Interpreter='latex') 
xlabel('t [T]',Interpreter='latex')
ylabel('$|a_{Cart}-a_{Gauss}|/a0 [-]$',Interpreter='latex')
xlim([0 K]) 
grid on; 
 
figure(8) %Eccentricity 
Ecc = tiledlayout(2,2);  
title(Ecc,'Eccentricity',Interpreter='latex') 
nexttile 
plot(t/T,y_kepGauss(:,2),'b-')
xlabel('t [T]',Interpreter='latex');ylabel('e [-]',Interpreter='latex')
hold on  
plot(t/T,y_kep(:,2),'r-')
xlabel('t [T]',Interpreter='latex');ylabel('e [-]',Interpreter='latex')
plot(t/T, E,'-','Color','[1, 0.84, 0]',LineWidth=2)
xlabel('t [T]',Interpreter='latex');ylabel('e [-]',Interpreter='latex')
legend('Gauss equations','Cartesian','Secular (filtered)') 
title([num2str(K) ' Periods'],Interpreter='latex') 
xlim([0 K]) 
grid on; 
hold off 
 
nexttile 
plot(t(1:n/10)/T,y_kepGauss(1:(n/10),2),'b-')
xlabel('t [T]',Interpreter='latex');ylabel('e [-]',Interpreter='latex')
hold on          
plot(t(1:n/10)/T,y_kep(1:(n/10),2),'r-')
xlabel('t [T]',Interpreter='latex');ylabel('e [-]',Interpreter='latex')
plot(t/T, E,'-','Color','[1, 0.84, 0]',LineWidth=2)
xlabel('t [T]',Interpreter='latex');ylabel('e [-]',Interpreter='latex')
legend('Gauss equations','Cartesian','Secular (filtered)') 
title([num2str(K/10) ' Periods'],Interpreter='latex') 
xlim([0 K/10]) 
grid on; 
hold off 
 
nexttile(3,[1 2]) 
semilogy(t/T,abs(y_kep(:,2)-y_kepGauss(:,2)),'b-') 
title('Relative error',Interpreter='latex')
xlabel('t [T]',Interpreter='latex')
ylabel('$|e_{Cart}-e_{Gauss}|/e0 [-]$',Interpreter='latex')
xlim([0 K]) 
grid on; 
 
figure(9) %Inclination 
Inc = tiledlayout(2,2);  
title(Inc,'Inclination',Interpreter='latex') 
nexttile 
plot(t/T,rad2deg(y_kepGauss(:,3)),'b-')
xlabel('t [T]',Interpreter='latex');ylabel('i [deg]',Interpreter='latex')
hold on  
plot(t/T,rad2deg(y_kep(:,3)),'r-')
xlabel('t [T]',Interpreter='latex');ylabel('i [deg]',Interpreter='latex')
plot(t/T, I,'-','Color','[1, 0.84, 0]',LineWidth=2)
xlabel('t [T]',Interpreter='latex');ylabel('i [deg]',Interpreter='latex')
title([num2str(K) ' Periods'],Interpreter='latex') 
legend('Gauss equations','Cartesian','Secular (filtered)') 
xlim([0 K]) 
grid on; 
hold off 
 
nexttile 
plot(t(1:n/10)/T,rad2deg(y_kepGauss(1:(n/10),3)),'b-')
xlabel('t [T]',Interpreter='latex');ylabel('i [deg]',Interpreter='latex')
hold on          
plot(t(1:n/10)/T,rad2deg(y_kep(1:(n/10),3)),'r-')
xlabel('t [T]',Interpreter='latex');ylabel('i [deg]',Interpreter='latex')
plot(t(1:n/10)/T, I(1:(n/10)),'-','Color','[1, 0.84, 0]',LineWidth=2)
xlabel('t [T]',Interpreter='latex');ylabel('i [deg]',Interpreter='latex')
legend('Gauss equations','Cartesian','Secular (filtered)') 
title([num2str(K/10) ' Periods'],Interpreter='latex') 
xlim([0 K/10]) 
grid on; 
hold off 
 
nexttile(3,[1 2]) 
semilogy(t/T,abs(y_kep(:,3)-y_kepGauss(:,3))/(2*pi),'b-') 
title('Relative error',Interpreter='latex') 
xlabel('t [T]',Interpreter='latex')
ylabel('$|i_{Cart}-i_{Gauss}|/2\pi [-]$',Interpreter='latex')
xlim([0 K]) 
grid on; 
 
figure(10) %OMEGA 
OM = tiledlayout(2,2); 
title(OM,'Right Ascension',Interpreter='latex') 
nexttile 
plot(t/T,rad2deg(y_kepGauss(:,4)),'b-')
xlabel('t [T]',Interpreter='latex');ylabel('$\Omega$ [deg]',Interpreter='latex')
hold on  
plot(t/T,rad2deg(unwrap(y_kep(:,4))),'r-')
xlabel('t [T]',Interpreter='latex');ylabel('$\Omega$ [deg]',Interpreter='latex')
plot(t/T, O,'-','Color','[1, 0.84, 0]',LineWidth=2)
xlabel('t [T]',Interpreter='latex');ylabel('$\Omega$ [deg]',Interpreter='latex')
title([num2str(K) ' Periods'],Interpreter='latex') 
legend('Gauss equations','Cartesian','Secular (filtered)') 
xlim([0 K]) 
grid on; 
hold off 
 
nexttile 
plot(t(1:n/10)/T,rad2deg(y_kepGauss(1:(n/10),4)),'b-')
xlabel('t [T]',Interpreter='latex');ylabel('$\Omega$ [deg]',Interpreter='latex')
hold on          
plot(t(1:n/10)/T,rad2deg(unwrap(y_kep(1:(n/10),4))),'r-')
xlabel('t [T]',Interpreter='latex');ylabel('$\Omega$ [deg]',Interpreter='latex')
plot(t(1:n/10)/T, O(1:n/10),'-','Color','[1, 0.84, 0]',LineWidth=2)
xlabel('t [T]',Interpreter='latex');ylabel('$\Omega$ [deg]',Interpreter='latex')
legend('Gauss equations','Cartesian','Secular (filtered)') 
title([num2str(K/10) ' Periods'],Interpreter='latex') 
xlim([0 K/10]) 
grid on; 
hold off 
 
nexttile(3,[1 2]) 
semilogy(t/T,abs(y_kep(:,4)-y_kepGauss(:,4))/(2*pi),'b-')
xlabel('t [T]',Interpreter='latex')
ylabel('$|\Omega_{Cart}-\Omega_{Gauss}|/2\pi [-]$',Interpreter='latex')
title('Relative error',Interpreter='latex') 
xlim([0 K]) 
grid on; 
 
% nexttile 
% plot(t/3600, y_kep(:,4),'b-',LineWidth=2) 
% hold on 
% plot(t/3600, y_aprrox(:,1),'r-',LineWidth=1) 
% legend('Gauss equations','Secular') 
% grid on; 
% hold off 
 
figure(11) %omega 
om = tiledlayout(2,2); 
title(om,'Argument of pericentre',Interpreter='latex') 
nexttile 
plot(t/T,rad2deg(y_kepGauss(:,5)),'b-')
xlabel('t [T]',Interpreter='latex');ylabel('$\omega$ [deg]',Interpreter='latex')
hold on  
plot(t/T,rad2deg(unwrap(y_kep(:,5))),'r-')
xlabel('t [T]',Interpreter='latex');ylabel('$\omega$ [deg]',Interpreter='latex')
plot(t/T, o,'-','Color','[1, 0.84, 0]',LineWidth=2)
xlabel('t [T]',Interpreter='latex');ylabel('$\omega$ [deg]',Interpreter='latex')
legend('Gauss equations','Cartesian','Secular (filtered)') 
title([num2str(K) ' Periods'],Interpreter='latex') 
xlim([0 K]) 
grid on; 
hold off 
 
nexttile 
plot(t(1:n/10)/T,rad2deg(y_kepGauss(1:(n/10),5)),'b-')
xlabel('t [T]',Interpreter='latex');ylabel('$\omega$ [deg]',Interpreter='latex')
hold on 
plot(t(1:n/10)/T,rad2deg(unwrap(y_kep(1:(n/10),5))),'r-')
xlabel('t [T]',Interpreter='latex');ylabel('$\omega$ [deg]',Interpreter='latex')
plot(t(1:n/10)/T, o(1:n/10),'-','Color','[1, 0.84, 0]',LineWidth=2)
xlabel('t [T]',Interpreter='latex');ylabel('$\omega$ [deg]',Interpreter='latex')
legend('Gauss equations','Cartesian','Secular (filtered)') 
title([num2str(K/10) ' Periods'],Interpreter='latex') 
xlim([0 K/10]) 
grid on; 
hold off

nexttile(3,[1 2])
semilogy(t/T,abs(y_kep(:,5)-y_kepGauss(:,5))/(2*pi),'b-')
xlabel('t [T]',Interpreter='latex')
ylabel('$|\omega_{Cart}-\omega_{Gauss}|/2\pi [-]$',Interpreter='latex')
title('Realtive error',Interpreter='latex')
xlim([0 K])
grid on;
% figure(17)
% plot(t/3600, o,'r-',LineWidth=2)% hold on
% plot(t/3600, y_aprrox(:,2),'b-',LineWidth=2)% grid on;
% hold off
figure(12) %theta
theta = tiledlayout(2,2);
title(theta,'True anomaly',Interpreter='latex')

nexttile
plot(t/T,rad2deg(y_kepGauss(:,6)),'b-')
xlabel('t [T]',Interpreter='latex');ylabel('f [deg]',Interpreter='latex')
hold on 
plot(t/T,rad2deg(unwrap(y_kep(:,6))),'r-')
xlabel('t [T]',Interpreter='latex');ylabel('f [deg]',Interpreter='latex')
plot(t/T, F,'-','Color','[1, 0.84, 0]',LineWidth=2)
xlabel('t [T]',Interpreter='latex');ylabel('f [deg]',Interpreter='latex')
legend('Gauss equations','Cartesian','Secular (filtered)',Interpreter='latex')
title([num2str(K) ' Periods'],Interpreter='latex')
xlim([0 K])
grid on;
hold off

nexttile
plot(t(1:n/10)/T,rad2deg(y_kepGauss(1:(n/10),6)),'b-')
xlabel('t [T]',Interpreter='latex');ylabel('f [deg]',Interpreter='latex')
hold on
plot(t(1:n/10)/T,rad2deg(unwrap(y_kep(1:(n/10),6))),'r-')
xlabel('t [T]',Interpreter='latex');ylabel('f [deg]',Interpreter='latex')
plot(t(1:n/10)/T, F(1:(n/10)),'-','Color','[1, 0.84, 0]',LineWidth=2)
xlabel('t [T]',Interpreter='latex');ylabel('f [deg]',Interpreter='latex')
legend('Gauss equations','Cartesian','Secular (filtered)')
title([num2str(K/10) ' Periods'],Interpreter='latex')
xlim([0 K/10])
grid on;
hold off

nexttile(3,[1 2])
semilogy(t/T,abs(unwrap(y_kep(:,6))-y_kepGauss(:,6))./abs(y_kepGauss(:,6)),'b-')
xlabel('t [T]',Interpreter='latex')
ylabel('$|f_{Cart}-f_{Gauss}|/f_{Gauss} [-]$',Interpreter='latex')
title('Relative error',Interpreter='latex')
xlim([0 K])
grid on;

%% TRACKING REAL OBJECT 43247 (debris)

% Modeling perturbations with Gauss equations
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
%angular velocity of Earth around its axis:
we = deg2rad(15.04)/3600; %[rad/s]
%planetary constant of the Earth:
mu_e = astroConstants(13); %[km^3/s^2]
%Radius of the Earth:
Re = astroConstants(23); %[km]
adebris = 15051; %semi-major axis [km]
edebris = 0.5449;  %eccentricity
idebris = deg2rad(54.9514); %inclination [rad]
% I have to choose the W,w,f0 parameteres:
OMEGA0debris = deg2rad(263.5589); %[rad]
omega0debris = deg2rad(139.9687); %[rad]
n=10000;
f0debris = deg2rad(0);   %[rad]
Tdebris = 2*pi*sqrt( adebris^3/mu_e); K = 60;
tspan = linspace(0,K*Tdebris,n);
   
% SRP + J2
J2 = astroConstants(9);
AU = astroConstants(2);
Cd = 2.1; % drag coefficient
AM = 0.0254; % spacecraft area-to-mass ratio relative to the air [m^2/kg]

% Propagation of keplerian elements of the debris
y0_kepGaussdebris = [adebris,edebris,idebris,OMEGA0debris,omega0debris,f0debris];
[r0debris, v0debris] = kep2cartmio(adebris,edebris,idebris,omega0debris,f0debris,OMEGA0debris,mu_e); y0_cartdebris = [r0debris;v0debris];
[tdebris,y_kepGaussdebris] = ode113(@(tdebris,ydebris) odefunperturb(tdebris,ydebris, mu_e, Re, J2,AM,Cd,we,'RSW') , tspan, y0_kepGaussdebris, options);
tempocomputGauss=toc;

%Let's filter the keplerian elements
Adebris = movmean(y_kepGaussdebris(:,1), 400);
Edebris = movmean(y_kepGaussdebris(:,2), 400);
Idebris = movmean(rad2deg(y_kepGaussdebris(:,3)), 400);
Odebris = movmean(rad2deg(y_kepGaussdebris(:,4)), 400);
odebris = movmean(rad2deg(y_kepGaussdebris(:,5)), 400);
Fdebris = movmean(rad2deg(unwrap(y_kepGaussdebris(:,6))), 400);

%10 days
tspaneph=linspace(0,K*Tdebris,11);
% [e,r_p,i,OMEGA,w,tp,n,M,theta,a,ra,tsid]
real_data=importdata('horizon2.mat');

figure()  %Semi-major axis   
plot(tdebris/Tdebris,y_kepGaussdebris(:,1),'r-')
hold on
plot(tdebris/Tdebris,Adebris,'-','Color','[1, 0.84, 0]',LineWidth=2)
plot(tspaneph/Tdebris,real_data(:,11),'k*')
title('Semi-major axis',[num2str(K) ' Periods'],Interpreter='latex')
xlabel('t [T]',Interpreter='latex');ylabel('a [km]',Interpreter='latex')
legend('Gauss','filtered','Ephemerides',Interpreter='latex')
grid on
hold off


figure() %Eccentricity   
plot(tdebris/Tdebris,y_kepGaussdebris(:,2),'b-')
hold on
plot(tdebris/Tdebris,Edebris,'-','Color','[1, 0.84, 0]',LineWidth=2)
plot(tspaneph/Tdebris,real_data(:,2),'k*')
title('Eccentricity debris',Interpreter='latex') 
xlabel('t [T]',Interpreter='latex');ylabel('e [-]',Interpreter='latex')
legend('Gauss','filtered','Ephemerides',Interpreter='latex')
grid on

figure() %Inclination   
plot(tdebris/Tdebris,rad2deg(y_kepGaussdebris(:,3)),'c-')
hold on
plot(tdebris/Tdebris,Idebris,'-','Color','[1, 0.84, 0]',LineWidth=2)
plot(tspaneph/Tdebris,real_data(:,4),'k*')
title('Inclination debris',Interpreter='latex')
xlabel('t [T]',Interpreter='latex');ylabel('i [deg]',Interpreter='latex')
legend('Gauss','filtered','Ephemerides',Interpreter='latex')
grid on

figure() %OMEGA  
plot(tdebris/Tdebris,rad2deg(y_kepGaussdebris(:,4)),'y-')
hold on
plot(tdebris/Tdebris,Odebris,'-','Color','[1, 0.84, 0]',LineWidth=2)
plot(tspaneph/Tdebris,real_data(:,5),'k*')
title('Right Ascension debris',Interpreter='latex') 
xlabel('t [T]',Interpreter='latex');ylabel('$\Omega$ [deg]',Interpreter='latex')
legend('Gauss','filtered','Ephemerides',Interpreter='latex')
grid on

figure() %omega  
plot(tdebris/Tdebris,rad2deg(y_kepGaussdebris(:,5)),'g-')
hold on
plot(tdebris/Tdebris,odebris,'-','Color','[1, 0.84, 0]',LineWidth=2)
plot(tspaneph/Tdebris,real_data(:,6),'k*')
title('Argument of pericentre debris',Interpreter='latex') 
xlabel('t [T]',Interpreter='latex');ylabel('$\omega$ [deg]',Interpreter='latex')
legend('Gauss','filtered','Ephemerides',Interpreter='latex')
grid on
