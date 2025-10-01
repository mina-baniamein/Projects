%function [n,p] = n_p_selector(sim_dat,t,theta0,Ts,model,p_cross)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation and Learning in Aerospace Project A.Y. 24-25 


% Outputs : n : [1 x 1] system order
%           p : [1x 1] past time window lenght

% Authors:  Alessandro Castelanelli Oddo (alessandro.castelanelli@polimi.it)
%            (@polimi.it)                     
%            (@polimi.it)                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Data
sim_dat = simulation_data;
t= (0:sample_time:simulation_time)';
p_cross = 0.75; 
Ts = sample_time;
theta0 = [Xu Xq Mu Mq Xd Md] + 0.5 * (-0.5 + rand(1,6)).*[Xu Xq Mu Mq Xd Md];
model_fun = @drone_model_PBSID;

% Data collection
u = sim_dat.Mtot; 
y = sim_dat.q;

% Vector definition
u_cross = u(ceil(length(u)*p_cross)+1:end);
u = u(1:ceil(length(u)*p_cross));
y_cross = y(ceil(length(y)*p_cross)+1:end);
y = y(1:ceil(length(y)*p_cross));
t_cross = t(ceil(length(t)*p_cross)+1:end);
t = t(1:ceil(length(t)*p_cross));

% Cutting frequency definition
% if vargin > 7
%     if varargin{1} == 1 
%         fc = 20;
%     else
%         fc = 10;
%     end
% else
%     fc = 15;
% end
fc = 20;
[u,y] = bettering_solution(u,y,fc,Ts);

% p vect definition
p_vect = [6,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26,28,30];
n_max = 4; i_max = 20;index = 1;
% PBSID cycle
error = zeros(n_max,i_max);
for i = 1 : i_max
    for n = 1 : n_max
        [D_PBSID,S_PBSID,V_PBSID,Y,N] = pbsid_1part(u,y,p_vect(i),1);
        [A_PBSID,B_PBSID,C_PBSID,K_PBSID] = pbsid_2part(D_PBSID,n,S_PBSID,V_PBSID,Y,N,p_vect(i),u,y);
        y_id = lsim(ss(A_PBSID,B_PBSID,C_PBSID,D_PBSID), u, t);
        [identification_PBSID,error_PBSID(index,:)] = greyest_structuring(u,y_id,Ts,theta0,model_fun,real_parameters);
        [A_PBSID,B_PBSID,C_PBSID,D_PBSID] = drone_model_grey(identification_PBSID.parameters,0);
        A = A_PBSID; B = B_PBSID; C([2,4],:) = C_PBSID; D([2,4],:) = D_PBSID;
        simulation_data_PBSID = sim('Simulator_Single_Axis','SrcWorkspace', 'current');
        y_id = simulation_data_PBSID.q;
        y_id_cross = y_id(ceil(length(y_id)*p_cross)+1:end);
        error(n,i) = max(abs((y_cross-y_id_cross)));
        index = index + 1;
    end
end

[n,i] = find(min(error));
p = p_vect(i);



















