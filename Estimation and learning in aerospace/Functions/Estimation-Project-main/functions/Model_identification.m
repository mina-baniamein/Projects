function [identification varargout] = Model_identification(u,y,guess,sample_time,model_fun,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation and Learning in Aerospace Project A.Y. 24-25 
% Function to identify by grey_box ID

% Authors:  Alessandro Castelanelli Oddo (alessandro.castelanelli@polimi.it)
%           Mina Baniamein (mina.baniamein@polimi.it)                     
%            (@polimi.it)                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Measured acceleration and pitch rate
%y = [data.q data.ax];
% Data ordering and bring those in frequency domain
sim_data = iddata(y, u, sample_time);
data_fd = fft(sim_data); % output of the simulation in the frequency domain

% Initial guess for the identification
sys_init = idgrey(model_fun, guess, 'c');

% Actual Model Identification
identification = struct;
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
