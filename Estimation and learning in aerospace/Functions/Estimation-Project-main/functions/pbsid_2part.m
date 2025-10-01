function [A,B,C,K] = pbsid_2part(D,n,S_svd,V_svd,Y,N,p,u,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation and Learning in Aerospace Project A.Y. 24-25 
% Second part of the PBSID for SISO after the n-order determination by SVD. 
% This code will evaluate the remaining matricies to identificate.

% Authors:  Alessandro Castelanelli Oddo (alessandro.castelanelli@polimi.it)
%            (@polimi.it)                     
%            (@polimi.it)                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% State estimation creation and input vector without past widow
X_est_n = sqrtm(S_svd(1:n, 1:n)) * V_svd(:, 1:n)';
u_C = u(p:N-1) ;
% C matrix estimation by LSQR
b_C = Y - D*u_C; % Right-hand side
A_C = X_est_n'; % Left-hand side
C = lsqr(A_C,b_C); % LSQR solution
% Innovation, right-hand and left-hand side initialization
e = zeros(N-2-p,1); % Innovation
b_ABK = zeros(2*(N-2-p),1); % Right-hand side
A_ABK = zeros(2*(N-2-p),n^2+2*n); % Left-hand side
A_a_mini = zeros(n,n^2); % Support matrix for A_ABK construction 
index = 1;
% Innovation, right-hand and left-hand side building
for i = p : N-2
    e(index) = y(i) - C' * X_est_n(:,i-p+1) - D* u(i);
    b_ABK (n*index-n+1 : n*index,1) = X_est_n(1:n,index+1);
    for j = 1 : n
        A_a_mini(j,j*n-n+1 : j*n) = X_est_n(:,index)';
        A_b_mini = eye(n)*u(i);
        A_k_mini = eye(n)*e(index);
    end
    A_ABK(n*index-n+1 : n*index, :) = [A_a_mini,A_b_mini,A_k_mini];
    index = index + 1;
 end
% A,B and K matricies estimation
ABK = lsqr(A_ABK, b_ABK);
% Solution reshaping
A = reshape(ABK(1:n^2),n,n);
B = ABK(n^2+1:n^2+n);
K = ABK(n^2+n+1:end);
A = A'; C = C' ;
end