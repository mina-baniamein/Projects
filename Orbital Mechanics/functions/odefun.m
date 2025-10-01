function dy = odefun(~,y,mu) 
% ode unperturbated 
% 
% PROTOTYPE:
%   dy=odefunperturb(t,y,mu,Re,J2,AM,Cd,we)
% 
% INPUT: 
%  t[1]       Time (can be omitted, as the system is autonomous)
%  y[6x1]     Cartesian state of the body (rx,ry,rz,vx,vy,vz) [L, L/T]
%  muP[1]     gravitational parameter of the primary [L^3/T^2]
% 
% OUTPUT:
%  dy[6x1]    Derivative of the state [L/T^2, L/T^3]
% 
% CONTRIBUTORS:
%  Pierpaolo Di Carlo 10767871
%  Alessandro Bellezza 10673485
%  Gaia Trovatelli 10582310
%  Mina Baniamein 10627453
%  
% VERSIONS
%  First versions

r=y(1:3);
v=y(4:6);
rnorm=norm(r);
dy=[v;(-mu/rnorm^3)*r];
end