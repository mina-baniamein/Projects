function [rcart,vcart] = kep2cartmio(a, ecc, inc, w, nu, RAAN,mu)
% this function returns r and v in cartesian coordinates starting forom keplerian elements
% 
% PROTOTYPE:
%   [rcart,vcart] = kep2cartmio(a, ecc, inc, w, nu, RAAN,mu)
% 
% INPUT: 
%  a[1]       semi major axis of the orbit [km]
%  ecc[1]     eccentricity 
%  inc[1]     inclination [rad]
%  w[1]       argment of perigee [rad]
%  nu[1]      RAAN (right ascension of the ascending node) [rad]
%  mu[1]      gravitational parameter of the primary [L^3/T^2]
% 
% OUTPUT:
%  [rcart,vcart]    cartesian vectors of radius and velocity
% 
% CONTRIBUTORS:
%  Pierpaolo Di Carlo 10767871
%  Alessandro Bellezza 10673485
%  Gaia Trovatelli 10582310
%  Mina Baniamein 10627453
%  
% VERSIONS
%  First version
        p = a*(1-ecc^2);
        h=sqrt(mu*p);
        r_0 = p/(1+ecc*cos(nu));
        x_perifocal=r_0*cos(nu);
        y_pericola=r_0*sin(nu);
% CALCOLO r nel perifocal        
        rperifocal=[x_perifocal,y_pericola,0]';
        Vx_=(mu/h)*(-sin(nu));
        Vy_=(mu/h)*(ecc+cos(nu));
% CALCOLO v nel perifocal        
        vperifocal=[Vx_,Vy_,0]';
% CREO LE MATRICI DI ROTAZIONE
        Rw=[ cos(w) sin(w) 0; -sin(w) cos(w) 0; 0 0 1 ];
        Ri=[ 1 0 0; 0 cos(inc) sin(inc); 0 -sin(inc) cos(inc) ];
        ROmega=[ cos(RAAN) sin(RAAN) 0; -sin(RAAN) cos(RAAN) 0; 0 0 1 ];
% OUTPUT:        
        rcart=(Rw*Ri*ROmega)'*rperifocal;
        vcart=(Rw*Ri*ROmega)'*vperifocal;
end
