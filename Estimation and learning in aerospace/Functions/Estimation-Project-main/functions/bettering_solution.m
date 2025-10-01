function [u_smooth,y_smooth] = bettering_solution(u_in,y_in,fc,Ts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation and Learning in Aerospace Project A.Y. 24-25 
% Function to reduce the noise, increasing the resolution of the solution

% Authors:  Alessandro Castelanelli Oddo (alessandro.castelanelli@polimi.it)
%           Mina Baniamein (mina.baniamein@polimi.it)                     
%            (@polimi.it)                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sampling frequency
fs = 1/Ts; 
[b,a] = butter(4, fc/(fs/2), 'low'); % Butterworth filter at 4 order
u_filt = filtfilt(b, a, u_in);
y_filt = filtfilt(b, a, y_in);

% Normalization (porta i dati tra 0 e 1)
u_norm = (u_filt - min(u_filt)) / (max(u_filt) - min(u_filt));
y_norm = (y_filt - min(y_filt)) / (max(y_filt) - min(y_filt));

% Frame length of Savitzky-Golay filter
framelen = fs/fc * 5; % 
framelen = round(framelen) + (mod(round(framelen), 2) == 0) * (2 * (round(framelen) < framelen) - 1);
framelen = 21;
% Smoothing with Savitzky-Golay filter
u_smooth = sgolayfilt(u_norm, 3, framelen);
y_smooth = sgolayfilt(y_norm, 3, framelen);

% Plotting
figure('Name','Optimizing Solution');
subplot(2,1,1);
plot(u_in, 'r--'); hold on; plot(u_smooth, 'b', 'LineWidth', 1.5);
legend('Original', 'Filtered and Normalized');
title('Input signal','FontSize',20);

subplot(2,1,2);
plot(y_in, 'r--'); hold on; plot(y_smooth, 'b', 'LineWidth', 1.5);
legend('Original', 'Filtered and Normalized');
title('Output signal','FontSize',20);

end