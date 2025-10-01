function [identification, varargout] = greyest_structuring(u,y,Ts,theta0,fun,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation and Learning in Aerospace Project A.Y. 24-25 
% Structuring maded by the use of the matlab function greyest

% Authors:  Alessandro Castelanelli Oddo (alessandro.castelanelli@polimi.it)
%            (@polimi.it)                     
%            (@polimi.it)                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data ordering and bring those in frequency domain
sim_data = iddata(y, u, Ts);
data_fd = fft(sim_data); % output of the simulation in the frequency domain

sys_init = idgrey(fun, theta0, 'c');
identification = struct;

% Greyest Options
options = greyestOptions('Display', 'on', 'SearchMethod', 'lsqnonlin');
options.SearchOptions.FunctionTolerance = 1e-6;
% Actual Model Identification
estimated_model = greyest(data_fd, sys_init);
identification.parameters = estimated_model.Report.Parameters.ParVector;
identification.fit = estimated_model.Report.Fit.FitPercent;
identification.fpe = estimated_model.Report.Fit.FPE;
identification.covariance = getcov(estimated_model);
identification.matrix={estimated_model.A; estimated_model.B; estimated_model.C; estimated_model.D};
identification.estimated_model=estimated_model;

% Data error
if nargout>=2
    real_parameters = varargin{1};
    varargout{1} = (identification.parameters-real_parameters) ./ real_parameters * 100; %Estimation Error
end
end