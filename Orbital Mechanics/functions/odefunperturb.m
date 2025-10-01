function dy = odefunperturb(~,y, mu, Re, J2,AM,Cd,we,frame)
% ode perturbated with J2 effect and air drag
% 
% PROTOTYPE:
%   dy=odefunperturb(t,y,mu,Re,J2,AM,Cd,we)
% 
% INPUT: 
%  t[1]       Time (can be omitted, as the system is autonomous)
%  y[6x1]     Cartesian state of the body (rx,ry,rz,vx,vy,vz) [L, L/T]
%  mu[1]      gravitational parameter of the primary [L^3/T^2]
%  Re[1]      Radius od the planet [L]
%  J2[1]      oblateness of the earth
%  AM[1]      ratio between the cross Area and the mass [L^2/KG]
%  Cd[1]      drag coefficient
%  we[1]      angular vecolity of the planet [rad/s]
%  frame(str)  the frame in which we are computing the perturbations:
%  'cart' or 'RSW'
% 
% OUTPUT:
%  dy[6x1]    Derivative of the state [L/T^2, L/T^3] for 'cart'
%  dy[6x1]    Derivative of the keplerian elements for 'RSW'
% 
% CONTRIBUTORS:
%  Pierpaolo Di Carlo 10767871
%  Alessandro Bellezza 10673485
%  Gaia Trovatelli 10582310
%  Mina Baniamein 10627453
%  
% VERSIONS
%  First versions
    switch frame
        case 'cart'
            r=y(1:3);
            v=y(4:6);
            rnorm=norm(r);
            h=rnorm-Re;
            wevector=[0;0;we];
            Vrel=v-cross(wevector,r);
            Vrelnorm=norm(Vrel);
            Vrelvers=Vrel/Vrelnorm;
            if h > 450 && h < 500
                h0=450;
                rho0=1.585e-12;
                H=60.828;
                rho=rho0*exp(-(h-h0)/H);
            elseif h > 500 && h < 600
                h0=500;
                rho0=6.967e-13;
                H=63.822;
                rho=rho0*exp(-(h-h0)/H);
            elseif h > 600 && h < 700
                h0=600;
                rho0=1.454e-13;
                H=71.835;
                rho=rho0*exp(-(h-h0)/H);
            elseif h > 700 && h < 800
                h0=700;
                rho0=3.614e-14;
                H=88.667;
                rho=rho0*exp(-(h-h0)/H);
            elseif h > 800 && h < 900
                h0=800;
                rho0=1.170e-14;
                H=124.64;
                rho=rho0*exp(-(h-h0)/H);
            elseif h > 900 && h < 1000
                h0=900;
                rho0=5.245e-15;
                H=181.05;
                rho=rho0*exp(-(h-h0)/H);
            elseif h > 1000
                h0=1000;
                rho0=3.019e-15;
                H=268;
                rho=rho0*exp(-(h-h0)/H);
            else
                rho=0;
            end
            %adrag=(-(1/2)*AM*Cd*rho*Vrelnorm^2*Vrelvers)*(10^3);
            adrag=(10^3)*(-(1/2)*AM*Cd*rho*(Vrelnorm^2)).*Vrelvers;
            %a_J2= (3/2)*(J2*mu*Re^2/rnorm^4)*((r(1)/rnorm)*(5*r(3)^2/rnorm^2-1) + (r(2)/rnorm)*(5*r(3)^2/rnorm^2-1) + (r(3)/rnorm)*(5*r(3)^2/rnorm^2-3));
            a_J2 = ((3/2)*(J2*mu*Re^2/rnorm^4)).*[(r(1)/rnorm)*(5*r(3)^2/rnorm^2 - 1); 
               (r(2)/rnorm)*(5*r(3)^2/rnorm^2 - 1); 
               (r(3)/rnorm)*(5*r(3)^2/rnorm^2 - 3)];
            dy= [v; (-mu/rnorm^3)*r+a_J2+adrag];
        case 'RSW'
            a = y(1); e = y(2); i = y(3); OMEGA = y(4); omega = y(5); f = y(6);
            b = a*sqrt(1-e^2); p = (b^2)/a; rr = p/(1+e*cos(f));h = sqrt(p*mu);
            [r,v]=kep2cartmio(a,e,i,omega,f,OMEGA,mu);
            rnorm=norm(r);
            altezza=rnorm-Re;
            wevector=[0;0;we];
            Vrel=v-cross(wevector,r);
            Vrelnorm=norm(Vrel);
            Vrelvers=Vrel/Vrelnorm;
            if altezza > 450 && altezza < 500
                h0=450;
                rho0=1.585e-12;
                H=60.828;
                rho=rho0*exp(-(altezza-h0)/H);
            elseif altezza > 500 && altezza < 600
                h0=500;
                rho0=6.967e-13;
                H=63.822;
                rho=rho0*exp(-(altezza-h0)/H);
            elseif altezza > 600 && altezza < 700
                h0=600;
                rho0=1.454e-13;
                H=71.835;
                rho=rho0*exp(-(altezza-h0)/H);
            elseif altezza > 700 && altezza < 800
                h0=700;
                rho0=3.614e-14;
                H=88.667;
                rho=rho0*exp(-(altezza-h0)/H);
            elseif altezza > 800 && altezza < 900
                h0=800;
                rho0=1.170e-14;
                H=124.64;
                rho=rho0*exp(-(altezza-h0)/H);
            elseif altezza > 900 && altezza < 1000
                h0=900;
                rho0=5.245e-15;
                H=181.05;
                rho=rho0*exp(-(altezza-h0)/H);
            elseif altezza > 1000
                h0=1000;
                rho0=3.019e-15;
                H=268;
                rho=rho0*exp(-(altezza-h0)/H);
            else
                rho=0;
            end
            %adrag=(-(1/2)*AM*Cd*rho*Vrelnorm^2*Vrelvers)*(10^3);
            adrag=(10^3)*(-(1/2)*AM*Cd*rho*(Vrelnorm^2)).*Vrelvers;
            %multiply with a rotation matrix
            Rotation=[cos(omega+f) sin(omega+f) 0;-sin(omega+f) cos(omega+f) 0;0 0 1]*[1 0 0;0 cos(i) sin(i);0 -sin(i) cos(i)]*[cos(OMEGA) sin(OMEGA) 0;-sin(OMEGA) cos(OMEGA) 0;0 0 1];
            adragRSW=Rotation*adrag;
            % Additional acceleration due to the J2
            a_J2 = -((3/2)*(J2*mu*Re^2/rr^4)).*[1-3*(sin(i)^2)*sin(f+omega)^2; 
               (sin(i)^2)*sin(2*(f+omega)); 
               sin(2*i)*sin(f + omega)];
            a_p=a_J2+adragRSW;
            % computing the derivatives 
            da = 2*((a^2)/h)*(e*sin(f)*a_p(1)+(p/rr)*a_p(2));
            de = (1/h)*(p*sin(f)*a_p(1)+((p+rr)*cos(f)+rr*e)*a_p(2));
            di = (rr*(cos(f+omega))/h)*a_p(3);
            dOMEGA = (rr*(sin(f+omega))/(h*sin(i))*a_p(3));
            domega = (1/(h*e))*(-p*cos(f)*a_p(1)+(p+rr)*sin(f)*a_p(2))-(rr*(sin(f+omega)*cos(i))/(h*sin(i))*a_p(3));
            df = (h/rr^2)+(1/(h*e))*(p*cos(f)*a_p(1)-(p+rr)*sin(f)*a_p(2));
            % Set the derivatives of the state
            dy = [da;de;di;dOMEGA;domega;df]; % faccio le derivate e le metto in un vettore [v; vDot]
    end
end








