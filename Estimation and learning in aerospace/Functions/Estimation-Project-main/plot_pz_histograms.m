function plot_pz_histograms(poles_grey, zeros_grey)
% plot_pz_histograms: Plots histograms of real and imaginary parts of poles and zeros
% with a dashed line at the mean value.
%
% Inputs:
%   poles_grey - complex vector of poles
%   zeros_grey - complex vector of zeros

figure('Name','Poles & Zeros Histograms','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

%% 1. Real part of poles
nexttile;
data = real(poles_grey);
histogram(data, 1000, 'FaceColor', [0.2 0.6 0.8]); % 30 bins, customizable color
hold on;
%xline(mean(data), 'r--', 'LineWidth', 2); % dashed line at mean
title('Poles - Real Part');
xlabel('Re'); ylabel('Count');
grid on;

%% 2. Imag part of poles
nexttile;
data = imag(poles_grey);
histogram(data, 1000, 'FaceColor', [0.2 0.6 0.8]);
hold on;
%xline(mean(data), 'r--', 'LineWidth', 2);
title('Poles - Imaginary Part');
xlabel('Im'); ylabel('Count');
grid on;

%% 3. Real part of zeros
nexttile;
data = real(zeros_grey);
histogram(data, 1000, 'FaceColor', [0.8 0.6 0.2]);
hold on;
%xline(mean(data), 'r--', 'LineWidth', 2);
title('Zeros - Real Part');
xlim([-1, 1]);
xlabel('Re'); ylabel('Count');
grid on;

%% 4. Imag part of zeros
nexttile;
data = imag(zeros_grey);
histogram(data, 1000, 'FaceColor', [0.8 0.6 0.2]);
hold on;
%xline(mean(data), 'r--', 'LineWidth', 2);
title('Zeros - Imaginary Part');
xlim([-2, 2]);
xlabel('Im'); ylabel('Count');
grid on;

end

