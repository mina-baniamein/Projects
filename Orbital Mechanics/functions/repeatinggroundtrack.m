function [T1,a1] = repeatinggroundtrack(mu,we,k,m)
% this function returns the period and the semi major orbit modified
% 
% PROTOTYPE:
%   [T1,a1]= repeatinggroundtrack(mu,we,k,m)
% 
% INPUT: 
%  muP[1]     gravitational parameter of the primary [L^3/T^2]
%  we[1]      earth angula velocity vector [rad/s]
%  k[1]       revolution of the satellite per orbit
%  m[1]       rotation of the earth per orbit
% 
% OUTPUT:
%  [T1,a1]    repeating period and modified semi major axis
% 
% CONTRIBUTORS:
%  Pierpaolo Di Carlo 10767871
%  Alessandro Bellezza 10673485
%  Gaia Trovatelli 10582310
%  Mina Baniamein 10627453
%  
% VERSIONS
%  First versions
T1=(2*pi*m)/(we*k);
a1=(mu/(we*k/m)^2)^(1/3);
end
