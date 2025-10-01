function [D,S_svd,V_svd,Y,N] = pbsid_1part(u,y,p,blocking_figure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation and Learning in Aerospace Project A.Y. 24-25 
% First part of the PBSID for SISO until the n-order determination by SVD. 
% This code generate an anti-diagonal matrix "Gamma_Delta_pp" with the 
% assumpion of the future window number is equal to the past (f == p)

% Input :   u : [Nx1] input signal
%           y : [Nx1] output signal
%           p : [1x1] past window
%           blocking_figure : generic parameter only to block the figure
%           generation

% Authors:  Alessandro Castelanelli Oddo (alessandro.castelanelli@polimi.it)
%            (@polimi.it)                     
%            (@polimi.it)                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the total number of data points
N = length(y);
% Output vector without past widow
Y = y(p:N-1);
Phi = zeros(N - p - 1, 2*p+1);  % Regression matrix inizialization
% Vector Z generation  
for k = p:N-1
    % Generate Z^(k-1, k-p) using u and y
    Z_kp_k = Z_k2_k1(u, y, k-1, k-p);  % Get p elements
    % Fill the regression matrix Phi
    Phi(k - p+1, 1:2*p) = Z_kp_k(:)';  % Flatten Z into a row
    Phi(k - p+1, end) = u(k);  % Last column is u(k)
end
% LSQR solution
x = lsqr(Phi,Y);
% Separation of results from LSQR
C_Delta_p = x(1:(length(x)-1)); 
D = x(end);
% Assumption on past and future window (f == p) and Hankel's matrix construction 
Gamma_Delta_pp = zeros(p,2*p);
C_Delta_p = [C_Delta_p', zeros(1,2*p)];
for i = 1:p
    Gamma_Delta_pp(i,:) = C_Delta_p(2*i-1: 2*i + 2*p -2);
end

% Full Z-array constuction for SVD
Z = [];
for k = 0 : N-p-1
    Z_temp = Z_k2_k1(u, y, p+k-1, k);
    Z = [Z, Z_temp];
end
% Singular Value Decomposition
[~, S_svd, V_svd] = svd(Gamma_Delta_pp * Z , 'econ');

% Useful for better_p funtion
if nargin < 4
    figure('Name','SVD Eig plottin')
    semilogy(diag(S_svd), 'o')
end

end
