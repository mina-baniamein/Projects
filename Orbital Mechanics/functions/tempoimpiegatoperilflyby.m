function [timeflyby]=tempoimpiegatoperilflyby(a1,a2,e1,e2,vinfmeno,vinfpiu,rsoi,rpericentro,mu)
% this function returns the time duration of the flyby
% 
% PROTOTYPE:
%   [timeflyby]=tempoimpiegatoperilflyby(a1,a2,e1,e2,vinfmeno,vinfpiu,rsoi,rpericentro,mu)
% 
% INPUT: 
%  a1[1]       semi major axis of the incoming hyperbola [km]
%  a2[1]       semi major axis of the outcoming hyperbola [km] 
%  e1[1]       ecccentricity if the incoming hyperbola [-]
%  e2[1]       eccentricity of the outcmoing hyperbola [-]
%  vinfmeno[1x3]    velocity vector of the incoming asymptote hyperbola 
%  vinfpiu[1x3]     velocity vector of the outcoming asymptote hyperbola 
%  rsoi[1]          radius of the sphere of influence of the planet [km]
%  rpericentro[1]   radius of the pericentre of the hyperbola[km]
%  mu[1]      gravitational parameter of the primary [L^3/T^2]
% 
% OUTPUT:
%  [timeflyby]    time duration of the flyby [DAYS] (considering a finite SOI)
% 
% CONTRIBUTORS:
%  Pierpaolo Di Carlo 10767871
%  Alessandro Bellezza 10673485
%  Gaia Trovatelli 10582310
%  Mina Baniamein 10627453
%  
% VERSIONS
%  First version
%E1=acosh((1/e2)*(1+rpericentro/abs(a2)));
E2=acosh((1/e2)*(1+rsoi/abs(a2)));
%E3=acosh((1/e1)*(1+rpericentro/abs(a1)));
E4=acosh((1/e1)*(1+rsoi/abs(a1)));
M1=e2*sinh(E2)-E2;
%M2=e2*sinh(E1)-E1;
%M3=e1*sinh(E3)-E3;
M4=e1*sinh(E4)-E4;
timeflyby=(a1^2*sqrt(e1^2-1)*(M4))/(rpericentro*sqrt(norm(vinfmeno)^2 + 2*mu/rpericentro))+(a2^2*sqrt(e2^2-1)*(M1))/(rpericentro*sqrt(norm(vinfpiu)^2+2*mu/rpericentro));
end