%% Assignment 1 group 2331

LAUNCHWINDOW=[2028,5,10,0,0,0; 2028,7,1,0,0,0];
GAWINDOW=[2028,7,31,0,0,0; 2028,9,1,0,0,0];
ARRIVALWINDOW=[2029,1,1,0,0,0; 2029,2,1,0,0,0];
au=astroConstants(2);
Asteroide=48;
n=1;
Rvenere=astroConstants(22);
altezzatmosfera=500; % [km]
guess=8000;
mu_v=astroConstants(12);
constraint=Rvenere+altezzatmosfera;
Planet=2; % flyby on Venus
StarterPlanet=1; % in this case we are launching from Mercurio
orbitType=0;
Nrev=0;
Ncase=0;
optionsLMR=0;

mjd2000lwi=date2mjd2000(LAUNCHWINDOW(1,1:6));
mjd2000lwf=date2mjd2000(LAUNCHWINDOW(2,1:6));

mjd2000gawi=date2mjd2000(GAWINDOW(1,1:6));
mjd2000gawf=date2mjd2000(GAWINDOW(2,1:6));

mjd2000awi=date2mjd2000(ARRIVALWINDOW(1,1:6));
mjd2000awf=date2mjd2000(ARRIVALWINDOW(2,1:6));

ascissa=(mjd2000lwi:n:mjd2000lwf);
terzadimensione=(mjd2000gawi:n:mjd2000gawf);
ordinata=(mjd2000awi:n:mjd2000awf);

dvtot=zeros(length(ascissa),length(terzadimensione),length(ordinata)); 
d=0;
tic
for mjd2000_d=mjd2000lwi:n:mjd2000lwf
    d=d+1;
    g=0;
    for mjd2000_ga=mjd2000gawi:n:mjd2000gawf
        g=g+1;
        h=0;
        if mjd2000_d > mjd2000_ga
            dvtot(d,g,:)=NaN;
            continue
        else
            [kepStarterPlanet,m_s] = uplanet(mjd2000_d,StarterPlanet); % kep = [a e i Om om theta] [km, rad]
            [riStarterPlanet,viStarterPlanet] = kep2cartmio(kepStarterPlanet(1),kepStarterPlanet(2),kepStarterPlanet(3),kepStarterPlanet(5),kepStarterPlanet(6),kepStarterPlanet(4),m_s); % (a, ecc, inc, w, nu, RAAN,mu)
            [kepPlanet,ks1] = uplanet(mjd2000_ga,Planet);
            [riPlanet,viPlanet] = kep2cartmio(kepPlanet(1),kepPlanet(2),kepPlanet(3),kepPlanet(5),kepPlanet(6),kepPlanet(4),ks1); % (a, ecc, inc, w, nu, RAAN,mu)
            TOF=(mjd2000_ga-mjd2000_d)*86400;
            [a_transfer,p_transfer,e_transfer,ERROR_transfer,Vi_transfer,Vf_transfer,TPAR_transfer,THETA_transfer] = lambertMR(riStarterPlanet,riPlanet,TOF,m_s,orbitType,Nrev,Ncase,optionsLMR);
            DELTA_Vdep=norm(Vi_transfer-viStarterPlanet');
            vinfmenooo=Vf_transfer-viPlanet';
        end
        for mjd2000_a=mjd2000awi:n:mjd2000awf
             h=h+1;
             if mjd2000_ga > mjd2000_a
                 dvtot(d,g,h)=NaN;
                 continue
             else
                 [kepAsteroide] = ephNEO(mjd2000_a,Asteroide); % kep = [a e i Om om theta] [km, rad]
                 TOF2=(mjd2000_a-mjd2000_ga)*86400;
                 [riAsteroide,viAsteroide] = kep2cartmio(kepAsteroide(1),kepAsteroide(2),kepAsteroide(3),kepAsteroide(5),kepAsteroide(6),kepAsteroide(4),m_s);
                 [a_2,p_2,e_2,ERROR_2,Vi_t2,Vf_t2,TPAR_transfer2,THETA_transfer2] = lambertMR(riPlanet,riAsteroide,TOF2,m_s,orbitType,Nrev,Ncase,optionsLMR);
                 vinfpiu=Vi_t2-viPlanet';
                 TURNANGLEE=acos((dot(vinfmenooo,vinfpiu)/(norm(vinfmenooo)*norm(vinfpiu))));
                 turnatot1=@(rpflyby) asin(1/(1+((rpflyby*(norm(vinfmenooo)^2))/mu_v)))+asin(1/(1+((rpflyby*(norm(vinfpiu)^2))/mu_v)))-TURNANGLEE;
                 rpflyby=fzero(turnatot1,guess);
                 if rpflyby <= constraint
                     dvtot(d,g,h)=NaN;
                     continue
                 else
                     Vphyperbolic1preflyby=sqrt(norm(vinfmenooo)^2+2*(mu_v/rpflyby));
                     Vphyperbolic2postflyby=sqrt(norm(vinfpiu)^2+2*(mu_v/rpflyby));
                     DELTA_VPga=norm(Vphyperbolic2postflyby-Vphyperbolic1preflyby);
                     DELTA_Varr=norm(Vf_t2-viAsteroide');
                     DELTA_VTOTALE= DELTA_VPga+DELTA_Varr+DELTA_Vdep;
                     dvtot(d,g,h)=DELTA_VTOTALE;
                 end
             end
        end
    end
end
computationaltimeofthe3FOR=toc;
mindV=min(min(min(dvtot)));

linearindex=find(dvtot==mindV);
[xmin,zmin,ymin]=ind2sub(size(dvtot),linearindex);

Timedeparture_min=ascissa(xmin);
Timeflyby_min=terzadimensione(zmin);
Timearrival_min=ordinata(ymin);
partenza=mjd20002date(Timedeparture_min);
dataflyby=mjd20002date(Timeflyby_min);
arrivo=mjd20002date(Timearrival_min);

timeofflight=(Timeflyby_min-Timedeparture_min)*86400;

[kepSP,ksunSP] = uplanet(Timedeparture_min,StarterPlanet); % departure from Mercury
[kepSPflyby,ksunSPflyby] = uplanet(Timeflyby_min,StarterPlanet);
[kepSPf,ksunSPf] = uplanet(Timearrival_min,StarterPlanet); 

[kepP,ksunP] = uplanet(Timeflyby_min,Planet); 
[kepPin,ksunPin] = uplanet(Timedeparture_min,Planet);
[kepPendmission,ksunPendmission] = uplanet(Timearrival_min,Planet); 

[kepAstin] = ephNEO(Timedeparture_min,Asteroide);
[kepAstfin] = ephNEO(Timearrival_min,Asteroide);
[kepAstflyby] = ephNEO(Timeflyby_min,Asteroide);

aSP=kepSP(1);
aP=kepP(1);
aASTEROIDE=kepAstin(1);

%Here we compute the position of the planets/asteroids in 3 differents
%time: departure, flyby and arrival

[riSP,viSP] = kep2cartmio(kepSP(1),kepSP(2),kepSP(3),kepSP(5),kepSP(6),kepSP(4),ksunSP);
[rPflyby,vPflyby] = kep2cartmio(kepP(1),kepP(2),kepP(3),kepP(5),kepP(6),kepP(4),ksunSP);
[rPendmission,vPendmission] = kep2cartmio(kepPendmission(1),kepPendmission(2),kepPendmission(3),kepPendmission(5),kepPendmission(6),kepPendmission(4),ksunSP);
[riP,viP] = kep2cartmio(kepPin(1),kepPin(2),kepPin(3),kepPin(5),kepPin(6),kepPin(4),ksunSP);
[rSPf,vSPf] = kep2cartmio(kepSPf(1),kepSPf(2),kepSPf(3),kepSPf(5),kepSPf(6),kepSPf(4),ksunSP);
[rSPflyby,vSPflyby] = kep2cartmio(kepSPflyby(1),kepSPflyby(2),kepSPflyby(3),kepSPflyby(5),kepSPflyby(6),kepSPflyby(4),ksunSP);
[rAstin,vAstin] = kep2cartmio(kepAstin(1),kepAstin(2),kepAstin(3),kepAstin(5),kepAstin(6),kepAstin(4),ksunSP);
[rAstfin,vAstfin] = kep2cartmio(kepAstfin(1),kepAstfin(2),kepAstfin(3),kepAstfin(5),kepAstfin(6),kepAstfin(4),ksunSP);
[rAstflyby,vAstflyby] = kep2cartmio(kepAstflyby(1),kepAstflyby(2),kepAstflyby(3),kepAstflyby(5),kepAstflyby(6),kepAstflyby(4),ksunSP);
TSP=2*pi*sqrt(aSP^3/ksunSP);
TP=2*pi*sqrt(aP^3/ksunSP);
TAst=2*pi*sqrt(aASTEROIDE^3/ksunSP);
y0SP=[riSP;viSP];
y0Pflyby=[rPflyby;vPflyby];
y0Pinitial=[riP;viP];
y0Astinitial=[rAstin;vAstin];
tspanSP=linspace(0,TSP,40000);
tspanP=linspace(0,TP,40000);
tspanAst=linspace(0,TAst,40000);

timeofflight2=(Timearrival_min-Timeflyby_min)*86400;

% Here we compute the 2 transfer arcs with the lambertMR function
[a_T1,p_T1,e_T1,ERROR_T1,Vi_T1,Vf_T1,TPAR_T1,THETA_T1] = lambertMR(riSP,rPflyby,timeofflight,ksunSP,orbitType,Nrev,Ncase,optionsLMR);
[a_T2,p_T2,e_T2,ERROR_T2,Vi_T2,Vf_T2,TPAR_T2,THETA_T2] = lambertMR(rPflyby,rAstfin,timeofflight2,ksunSP,orbitType,Nrev,Ncase,optionsLMR);
[a_transfer1,e_transfer1,i_transfer1,Omega_transfer1,omega_transfer1,theta_transfer1] = cart2kep(riSP, Vi_T1', ksunSP, 'deg');
[a_transfer2,e_transfer2,i_transfer2,Omega_transfer2,omega_transfer2,theta_transfer2] = cart2kep(rPflyby, Vi_T2', ksunSP, 'deg');
yTRANSFER10=[riSP;Vi_T1'];
yTRANSFER20=[rPflyby;Vi_T2'];
tspanARC1=linspace(0,timeofflight,40000);
tspanARC2=linspace(0,timeofflight2,40000);
tspantotal=linspace(0,timeofflight2+timeofflight,40000);
timeofthemission=(timeofflight+timeofflight2)/86400;

vInfplus=Vi_T2-vPflyby';
vInfminus=Vf_T1-vPflyby';

deltavfb=norm(vInfplus-vInfminus); % delta velocity given by the entire flyby
Turnangle=acos((dot(vInfminus,vInfplus)/(norm(vInfminus)*norm(vInfplus))));
Turnangletotale1=@(rpericentre) asin(1/(1+((rpericentre*(norm(vInfminus)^2))/mu_v)))+asin(1/(1+((rpericentre*(norm(vInfplus)^2))/mu_v)))-Turnangle;
rpericentre=fzero(Turnangletotale1,constraint);% this variable 'constraint' takes into account the radius of Venus and its atmosphere

Vphyperbolic1=sqrt(norm(vInfminus)^2+2*(mu_v/rpericentre)); % hyperbolic velocity for the incoming hyperbola at the pericentre
Vphyperbolic2=sqrt(norm(vInfplus)^2+2*(mu_v/rpericentre)); % hyperbolic velocity for the incoming hyperbola at the pericentre
deltavelocityga=norm(Vphyperbolic2-Vphyperbolic1);% delta velocity that I need to provid in order to combine the incoming and the outcoming hyperbola
deltavelocitydepaprture=norm(Vi_T1-viSP');% delta velocity needed for the departure from Mercury
deltvelocityarrive=norm(Vf_T2-vAstfin');% delta velocity needed for the arrival at the Asteroid
a1=-mu_v/norm(vInfminus)^2; % semi major axis of the incoming hyperbola
a2=-mu_v/norm(vInfplus)^2; % semi major axis of the outcoming hyperbola
delta1=asin(1/(1+((rpericentre*(norm(vInfminus)^2))/mu_v)));
delta2=asin(1/(1+((rpericentre*(norm(vInfplus)^2))/mu_v)));
e1=1/sin(delta1/2);
e2=1/sin(delta2/2);
b1=abs(a1)*sqrt(1-e1^2);
b2=abs(a2)*sqrt(1-e2^2);
rsoi=(0.72*au)*((mu_v/ksunP)^(2/5));

altitudeoftheclosestapproach=rpericentre-Rvenere; % altitude of the closest approach > atmosphere of Venus
[tempodelflyby]=tempoimpiegatoperilflyby(a1,a2,e1,e2,vInfminus,vInfplus,rsoi,rpericentre,mu_v);
timeperflyby=tempodelflyby/86400; % time duration of the flyby

rpflybyettore=[rpericentre;0;0];
Vhyperbolic1vettore=[0;Vphyperbolic1;0];
Vhyperbolic2vettore=[0;Vphyperbolic2;0];

yplanetoc1=[rpflybyettore;Vhyperbolic1vettore];
yplanetoc2=[rpflybyettore;-Vhyperbolic2vettore];

tspan=linspace(0,2000,1000);
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[tbplanet, ybplanet]=ode113(@(t,y) odefun(t,y,mu_v),tspan,yplanetoc1,options);
[taplanet, yaplanet]=ode113(@(t,y) odefun(t,y,mu_v),tspan,yplanetoc2,options);

% ORBIT PROPAGATION
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[tARC1,yARC1]=ode113(@(t,y) odefun(t,y,ksunSP),tspanARC1,yTRANSFER10,options); % transfer 1 
[tARC2,yARC2]=ode113(@(t,y) odefun(t,y,ksunSP),tspanARC2,yTRANSFER20,options); % transfer 2
[tStarterPlanet, yStarterPlanet]=ode113(@(t,y) odefun(t,y,ksunSP),tspanSP,y0SP,options); % Mercury orbit
[tPlanet, yPlanet]=ode113(@(t,y) odefun(t,y,ksunSP),tspanP,y0Pflyby,options); % Venus orbit
[tAst, yAst]=ode113(@(t,y) odefun(t,y,ksunSP),tspanAst,y0Astinitial,options); % Asteroid orbit
[tARCSP,yARCSP]=ode113(@(t,y) odefun(t,y,ksunSP),tspantotal,y0SP,options); % Mercury arc from the departure to the arrival (more than 1 orbit around the Sun)
[tARCP,yARCP]=ode113(@(t,y) odefun(t,y,ksunSP),tspantotal,y0Pinitial,options); % Venus arc from the departure to the arrival (more than 1 orbit around the Sun)
[tARCAst,yARCAst]=ode113(@(t,y) odefun(t,y,ksunSP),tspantotal,y0Astinitial,options); % Asteroid arc form the departure to arrival

figure()
plot3(yStarterPlanet(:,1),yStarterPlanet(:,2),yStarterPlanet(:,3),'r','LineStyle','--','Linewidth',3) % Mercury orbit
hold on
plot3(yARC1(:,1),yARC1(:,2),yARC1(:,3),'cyan','LineStyle','-','Linewidth',3) % transfer 1
plot3(yARC2(:,1),yARC2(:,2),yARC2(:,3),'cyan','LineStyle','-','Linewidth',3) % transfer 2
plot3(yPlanet(:,1),yPlanet(:,2),yPlanet(:,3),'b','LineStyle','--','Linewidth',3) % Venus orbit
plot3(yAst(:,1),yAst(:,2),yAst(:,3),'g','LineStyle','--','Linewidth',3) % Asteroid orbit
plot3(yARCAst(:,1),yARCAst(:,2),yARCAst(:,3),'g','LineStyle','-','Linewidth',3) % Asteroid arc form the departure to arrival
plot3(yARCSP(:,1),yARCSP(:,2),yARCSP(:,3),'r','LineStyle','-','Linewidth',3) % Mercury arc from the departure to arrival
plot3(yARCP(:,1),yARCP(:,2),yARCP(:,3),'b','LineStyle','-','Linewidth',3) % Venus arc from the departure to arrival
[Mercury] = plotPlanet(1, [riSP(1),riSP(2),riSP(3)], gca, 10); % Mercury position at the departure
[Venus] = plotPlanet(2, [riP(1),riP(2),riP(3)], gca, 10); % Venus position at the departure
scatter3(rAstin(1),rAstin(2),rAstin(3),60,'black','filled') % Asteroid position at the departure
[Venus] = plotPlanet(2, [rPendmission(1),rPendmission(2),rPendmission(3)], gca, 10); % Venus position at the end of the mission
[Venus] = plotPlanet(2, [rPflyby(1),rPflyby(2),rPflyby(3)], gca, 10); % Venus position at flyby
scatter3(rAstfin(1),rAstfin(2),rAstfin(3),60,'black','filled') % Asteroid position at the end of the mission
scatter3(rAstflyby(1),rAstflyby(2),rAstflyby(3),60,'black','filled') % Asteroid position at the flyby
[Sun] = plotPlanet(10, [0 0 0], gca, 17);
[Mercury] = plotPlanet(1, [rSPf(1),rSPf(2),rSPf(3)], gca, 10); % Mercury position at the end of the mission
[Mercury] = plotPlanet(1, [rSPflyby(1),rSPflyby(2),rSPflyby(3)], gca, 10); % Mercury position at the flyby
xlabel('x [Km]',Interpreter='latex');ylabel('y [km]',Interpreter='latex');zlabel('z [km]',Interpreter='latex');
title('Interplanetary mission',Interpreter='latex')
legend('Mercury orbit','Transfer arc 1','Transfer arc 2','Venus orbit','Asteroid orbit','Asteroid orbit arc during mission','Mercury motion during mission','Venus motion during mission','Asteroid',Interpreter='latex')
axis equal;
grid on
hold off

figure()
planet3D('Venus',struct('Units','km'));
hold on
plot3(ybplanet(:,1),ybplanet(:,2),ybplanet(:,3),'red','Linewidth',3)
plot3(yaplanet(:,1),yaplanet(:,2),yaplanet(:,3),'blue','Linewidth',3)
xlabel('x [Km]',Interpreter='latex');ylabel('y [km]',Interpreter='latex');zlabel('z [km]',Interpreter='latex');
title('Fly by hyperbola',Interpreter='latex')
legend('Venus','Incoming hyperbola','Outcoming hyperbola',Interpreter='latex')
axis equal;
grid on
hold off

ratioo=deltavelocityga/deltavfb;
fprintf('\n');
% VALUES
fprintf('Computational time for the 3 FOR cycles: %.7f sec \n',computationaltimeofthe3FOR);
fprintf('DeltaVtot_minimum= %.3f km/s \n',mindV);
fprintf('Departure day: %s \n',mat2str(partenza));
fprintf('Date of the flyby %s \n',mat2str(dataflyby));
fprintf('Arrival date: %s \n',mat2str(arrivo));
fprintf('DeltaV_ga: %.3f km/s \n',deltavelocityga);
fprintf('DeltaV_departure: %.3f km/s \n',deltavelocitydepaprture);
fprintf('DeltaV_arrival: %.3f km/s \n',deltvelocityarrive);
fprintf('Altitude of the closest approach: %.3f km \n',altitudeoftheclosestapproach);
fprintf('Flyby pericentre radius: %.2f km \n',altitudeoftheclosestapproach+Rvenere);
fprintf('Time duration of the flyby (considering a finite SOI): %.3f days \n',timeperflyby);
fprintf('Total time of flight: %.3f days \n',timeofthemission);
fprintf('Total time of flight (hours): %.3f h \n',timeofthemission*24);
fprintf('Ratio between the DV of the impulsive manoeuvre at the pericentre of the flyby and the total DV provided by the flyby: %.3f \n',ratioo);


