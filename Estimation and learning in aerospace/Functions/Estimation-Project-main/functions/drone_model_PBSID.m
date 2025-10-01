function [A,B,C,D] = drone_model_PBSID(params,Ts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation and Learning in Aerospace Project A.Y. 24-25 
% Function to describe the drone model

% Authors:  Alessandro Castelanelli Oddo (alessandro.castelanelli@polimi.it)
%            (@polimi.it)                     
%            (@polimi.it)                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters definition
Xu = params(1);
Xq = params(2);
Mu = params(3);
Mq = params(4);
Xd = params(5);
Md = params(6);

% Model definition
A=[Xu, Xq, -9.81; Mu, Mq, 0; 0, 1, 0];
B=[Xd; Md; 0];
C=[0, 1, 0]; 
D=0;

end