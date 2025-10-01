function [alpha, delta, lon, lat] = groundTrack(t0,tspan,y0,theta0,we,mu)
% this function returns the riht ascension, the declination, the longitude and the latitude
% of the satellite
% 
% PROTOTYPE:
%   [alpha, delta, lon, lat] = groundTrack(t0,tspan,y0,theta0,we,mu)
% 
% INPUT: 
%  t0[1]           initial time [s]
%  tspan[1xN]      time vector [s]
%  y0[3x6]         initial conditions [km,km\s]
%  theta0[1]       initial position of Greenwhich [rad]
%  we[1]           earth angula velocity vector [rad/s]
%  mu[1]           gravitational parameter of the primary [L^3/T^2]
% 
% OUTPUT:
%  [alpha, delta, lon, lat]      riht ascension, the declination, the longitude and the latitude of the satellite
% 
% CONTRIBUTORS:
%  Pierpaolo Di Carlo 10767871
%  
% VERSIONS
%  First versions
        options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
        [t, y]=ode113(@(t,y) odefun(t,y,mu),tspan,y0,options);
        r=y(:,1:3);
        R=vecnorm(r,2,2);
        theta=theta0 + we*(t-t0);
        alpha=atan2(r(:,2),r(:,1));
        delta=asin(r(:,3)./R);
        lat=wrapTo180(rad2deg(delta));
        lon=wrapTo180(rad2deg(alpha-theta));
end