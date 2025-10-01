function [mc_results, mc_summary] = monte_carlo_greybox( noise, input_mc, ...
    Xu, Xq, Mu, Mq, Xd, Md, sample_time, varargin)
%MONTE_CARLO_GREYBOX_IDENTIFICATION Performs Monte Carlo simulation for grey-box parameter identification
%
% Syntax:
%   [mc_results, mc_summary] = monte_carlo_greybox_identification(u_3ord, y_q, y_ax, ...
%       Xu, Xq, Mu, Mq, Xd, Md, sample_time)
%   [mc_results, mc_summary] = monte_carlo_greybox_identification(..., 'Name', Value)
%
% Inputs:
%   u_3ord      - Control input vector
%   y_q         - Pitch rate measurements
%   y_ax        - Longitudinal acceleration measurements
%   Xu,Xq,Mu,Mq,Xd,Md - True parameter values for reference
%   sample_time - Sampling time [s]
%
% Optional Name-Value Pairs:
%   'num_monte_carlo' - Number of Monte Carlo runs (default: 100)
%   'noise_enabler'   - Enable/disable noise (default: 1)
%   'pos_noise_std'   - Position noise std dev [m] (default: 0.0011)
%   'vel_noise_std'   - Velocity noise std dev [m/s] (default: 0.01)
%   'att_noise_std'   - Attitude noise std dev [rad] (default: deg2rad(0.33))
%   'rate_noise_std'  - Angular rate noise std dev [rad/s] (default: deg2rad(1))
%   'model_function'  - Model function handle (default: @drone_model_grey)
%   'verbose'         - Display progress (default: true)
%   'save_results'    - Save results to file (default: true)
%   'plot_results'    - Generate plots (default: true)
%
% Outputs:
%   mc_results  - Structure containing all Monte Carlo results
%   mc_summary  - Structure containing summary statistics
%
% Example:
%   [results, summary] = monte_carlo_greybox_identification(u_input, pitch_rate, ...
%       accel_x, -0.5, 2.1, -1.2, -8.5, 0.2, -0.8, 0.01, 'num_monte_carlo', 50);

%% Input Parsing
p = inputParser;

addRequired(p, 'Xu', @isnumeric);
addRequired(p, 'Xq', @isnumeric);
addRequired(p, 'Mu', @isnumeric);
addRequired(p, 'Mq', @isnumeric);
addRequired(p, 'Xd', @isnumeric);
addRequired(p, 'Md', @isnumeric);
addRequired(p, 'sample_time', @(x) isnumeric(x) && x > 0);

addParameter(p, 'num_monte_carlo', 100, @(x) isnumeric(x) && x > 0);
addParameter(p, 'noise_enabler', 1, @(x) isnumeric(x) && (x == 0 || x == 1));
addParameter(p, 'pos_noise_std', 0.0011, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'vel_noise_std', 0.01, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'att_noise_std', deg2rad(0.33), @(x) isnumeric(x) && x >= 0);
addParameter(p, 'rate_noise_std', deg2rad(1), @(x) isnumeric(x) && x >= 0);
addParameter(p, 'model_function', @drone_model_grey, @(x) isa(x, 'function_handle'));
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'save_results', true, @islogical);
addParameter(p, 'plot_results', true, @islogical);

parse(p, Xu, Xq, Mu, Mq, Xd, Md, sample_time, varargin{:});

% Extract parsed values
num_monte_carlo = p.Results.num_monte_carlo;
model_fun = p.Results.model_function;
verbose = p.Results.verbose;
save_results = p.Results.save_results;
plot_results = p.Results.plot_results;

%% Noise Configuration (as per existing definition)

noise.Enabler = 1;

noise.pos_stand_dev = noise.Enabler * 0.0011;                            	%[m]

noise.vel_stand_dev = noise.Enabler * 0.01;                               %[m/s]

noise.attitude_stand_dev = noise.Enabler * deg2rad(0.33);                 %[rad]
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(1);                   %[rad/s]

% Delays

delay.position_filter = 1;
delay.attitude_filter = 1;
delay.mixer = 1;

%% overwrite noise parameters
noise.Enabler = p.Results.noise_enabler;
noise.pos_stand_dev = noise.Enabler * p.Results.pos_noise_std;        % [m]
noise.vel_stand_dev = noise.Enabler * p.Results.vel_noise_std;        % [m/s]  
noise.attitude_stand_dev = noise.Enabler * p.Results.att_noise_std;   % [rad]
noise.ang_rate_stand_dev = noise.Enabler * p.Results.rate_noise_std;  % [rad/s]

%% Initialize Results Storage
mc_results = struct();
mc_results.identified_params = zeros(num_monte_carlo, 6);  % [Xu Xq Mu Mq Xd Md]
mc_results.errors = zeros(num_monte_carlo, 6);
mc_results.A_matrices = cell(num_monte_carlo, 1);
mc_results.B_matrices = cell(num_monte_carlo, 1);
mc_results.C_matrices = cell(num_monte_carlo, 1);
mc_results.D_matrices = cell(num_monte_carlo, 1);
mc_results.convergence_flags = zeros(num_monte_carlo, 1);
mc_results.noise_realizations = struct();

%% Display Configuration
if verbose
    fprintf('Starting Monte Carlo Grey-Box Parameter Identification...\n');
    fprintf('Noise Configuration:\n');
    fprintf('  Noise Enabler: %d\n', noise.Enabler);
    fprintf('  Position std dev: %.6f m\n', noise.pos_stand_dev);
    fprintf('  Velocity std dev: %.4f m/s\n', noise.vel_stand_dev);
    fprintf('  Attitude std dev: %.4f rad (%.2f deg)\n', noise.attitude_stand_dev, rad2deg(noise.attitude_stand_dev));
    fprintf('  Angular rate std dev: %.4f rad/s (%.2f deg/s)\n', noise.ang_rate_stand_dev, rad2deg(noise.ang_rate_stand_dev));
    fprintf('  Monte Carlo runs: %d\n\n', num_monte_carlo);
end

%% Real Parameter Vector (reference)
real_parameters = [Xu; Xq; Mu; Mq; Xd; Md];
ExcitationM = input_mc;
t=ExcitationM(:,1);

A=[Xu, Xq, -9.81; Mu, Mq, 0; 0, 1, 0];

B=[Xd; Md; 0];

C=[1, 0, 0; 0, 1, 0; 0, 0, 1; Xu, Xq, 0]; 

D=[0; 0; 0; Xd];

simulation_time=t(end)-t(1);

parameters_controller 

%% Monte Carlo Loop
if verbose
    tic;  % Start timing
end

for mc_run = 1:num_monte_carlo
    
    if verbose && mod(mc_run, 10) == 0
        fprintf('Progress: %d/%d runs completed\n', mc_run, num_monte_carlo);
    end
    
        % --- Generate random seeds for each noise block ---
    NoiseSeed_pos   = randi(1e6);
    NoiseSeed_vel   = randi(1e6);
    NoiseSeed_theta = randi(1e6);
    NoiseSeed_q     = randi(1e6);

    % Assign them to the base workspace (so Simulink can see them)
    assignin('base','NoiseSeed_pos',NoiseSeed_pos);
    assignin('base','NoiseSeed_vel',NoiseSeed_vel);
    assignin('base','NoiseSeed_theta',NoiseSeed_theta);
    assignin('base','NoiseSeed_q',NoiseSeed_q);

    % --- Run simulation ---
    simulation_data = sim('Simulator_Single_Axis_MC','SrcWorkspace','current');

    % Extract & identify just like above
    u_3ord = simulation_data.Mtot; 
    y_q    = simulation_data.q;
    y_ax   = simulation_data.ax;
   

    y_grey = [y_q y_ax];
    model_fun  = @drone_model_grey;

    %% Parameter Identification with Noisy Data
    % Generate new random guess for each run
    guess = [Xu Xq Mu Mq Xd Md] + 0.1 * (-0.5 + rand(1,6)).*[Xu Xq Mu Mq Xd Md];
    
    % Function call with noisy data
    try
        [identification_grey, error_grey] = Model_identification(u_3ord, y_grey, ...
                                                               guess, sample_time, model_fun, real_parameters);
        
        % Extract results
        A_grey = identification_grey.matrix{1};
        B_grey = identification_grey.matrix{2};
        C_grey = identification_grey.matrix{3};
        D_grey = identification_grey.matrix{4};
        
        % Store results - convert column vector to row vector
        mc_results.identified_params(mc_run, :) = identification_grey.parameters';
        mc_results.errors(mc_run, :) = error_grey';
        mc_results.A_matrices{mc_run} = A_grey;
        mc_results.B_matrices{mc_run} = B_grey;
        mc_results.C_matrices{mc_run} = C_grey;
        mc_results.D_matrices{mc_run} = D_grey;
        mc_results.convergence_flags(mc_run) = 1;  % Successful identification
        
    catch ME
        % Handle identification failures
        if verbose
            fprintf('Warning: Run %d failed with error: %s\n', mc_run, ME.message);
        end
        mc_results.identified_params(mc_run, :) = NaN(1, 6);
        mc_results.errors(mc_run) = NaN;
        mc_results.A_matrices{mc_run} = [];
        mc_results.B_matrices{mc_run} = [];
        mc_results.C_matrices{mc_run} = [];
        mc_results.D_matrices{mc_run} = [];
        mc_results.convergence_flags(mc_run) = 0;  % Failed identification
    end
    

end

if verbose
    elapsed_time = toc;
    fprintf('\nMonte Carlo simulation completed in %.2f seconds\n', elapsed_time);
else
    elapsed_time = 0;
end

%% Statistical Analysis
successful_runs = find(mc_results.convergence_flags == 1);
num_successful = length(successful_runs);
success_rate = num_successful / num_monte_carlo * 100;

if verbose
    fprintf('\nMonte Carlo Results Summary:\n');
    fprintf('Success rate: %.1f%% (%d/%d runs)\n', success_rate, num_successful, num_monte_carlo);
end

%% Initialize summary structure
mc_summary = struct();
mc_summary.num_runs = num_monte_carlo;
mc_summary.noise_config = noise;
mc_summary.success_rate = success_rate;
mc_summary.computation_time = elapsed_time;

if num_successful > 0
    % Calculate statistics only for successful runs
    successful_params = mc_results.identified_params(successful_runs, :);
    successful_errors = mc_results.errors(successful_runs);
    
    % Parameter statistics
    param_mean = mean(successful_params, 1);
    param_std = std(successful_params, 0, 1);
    param_names = {'Xu', 'Xq', 'Mu', 'Mq', 'Xd', 'Md'};
    
    if verbose
        fprintf('\nParameter Statistics (successful runs only):\n');
        fprintf('Parameter | True Value | Mean      | Std Dev   | Bias      | CV(%%)\n');
        fprintf('----------|------------|-----------|-----------|-----------|-------\n');
        
        for i = 1:6
            true_val = real_parameters(i);
            bias = param_mean(i) - true_val;
            cv = abs(param_std(i) / param_mean(i)) * 100;  % Coefficient of variation
            
            fprintf('%-9s | %10.4f | %9.4f | %9.4f | %9.4f | %6.2f\n', ...
                param_names{i}, true_val, param_mean(i), param_std(i), bias, cv);
        end
    end
    
    % Error statistics
    error_mean = mean(successful_errors);
    error_std = std(successful_errors);
    error_min = min(successful_errors);
    error_max = max(successful_errors);
    
    if verbose
        fprintf('\nIdentification Error Statistics:\n');
        fprintf('Mean error: %.6f\n', error_mean);
        fprintf('Std deviation: %.6f\n', error_std);
        fprintf('Min error: %.6f\n', error_min);
        fprintf('Max error: %.6f\n', error_max);
        
        %% Noise Impact Analysis
        fprintf('\nNoise Impact Analysis:\n');
        fprintf('Angular rate noise impact on pitch rate (y_q):\n');
        fprintf('  SNR: %.2f dB\n', 20*log10(rms(y_q)/noise.ang_rate_stand_dev));
        fprintf('Acceleration noise impact on longitudinal accel (y_ax):\n');
        fprintf('  Equivalent SNR: %.2f dB\n', 20*log10(rms(y_ax)/(noise.vel_stand_dev/sample_time)));
    end
    
    % Store statistics in summary
    mc_summary.parameter_statistics.mean = param_mean;
    mc_summary.parameter_statistics.std = param_std;
    mc_summary.parameter_statistics.true_values = real_parameters';
    mc_summary.error_statistics.mean = error_mean;
    mc_summary.error_statistics.std = error_std;
    mc_summary.error_statistics.min = error_min;
    mc_summary.error_statistics.max = error_max;
    
    %% Visualization
    if plot_results
        % Parameter distribution plots
        figure('Name', 'Monte Carlo Parameter Distributions', 'Position', [100 100 1200 800]);
        
        for i = 1:6
            subplot(2, 3, i);
            histogram(successful_params(:, i), 20, 'Normalization', 'probability');
            hold on;
            % Plot true value line
            ylim_vals = ylim;
            plot([real_parameters(i) real_parameters(i)], ylim_vals, 'r--', 'LineWidth', 2);
            % Plot mean line
            plot([param_mean(i) param_mean(i)], ylim_vals, 'g-', 'LineWidth', 2);
            
            xlabel(param_names{i});
            ylabel('Probability');
            title(sprintf('%s Distribution\nTrue=%.3f, Mean=%.3f±%.3f', ...
                param_names{i}, real_parameters(i), param_mean(i), param_std(i)));
            legend('Histogram', 'True Value', 'Mean', 'Location', 'best');
            grid on;
        end
        
        % Error distribution and noise correlation analysis
        figure('Name', 'Monte Carlo Error and Noise Analysis', 'Position', [200 200 1000 800]);
        
        subplot(2, 2, 1);
        histogram(successful_errors, 30);
        xlabel('Identification Error');
        ylabel('Frequency');
        title(sprintf('Error Distribution (Mean: %.4f ± %.4f)', error_mean, error_std));
        grid on;
        
        subplot(2, 2, 2);
        plot(successful_runs, successful_errors, 'b.-', 'MarkerSize', 8);
        xlabel('Monte Carlo Run Number');
        ylabel('Identification Error');
        title('Error vs Run Number');
        grid on;
        
        % % Noise realization examples
        % subplot(2, 2, 3);
        % run_example = min(10, num_successful);  % Show first 10 runs or available
        % for i = 1:run_example
        %     plot(mc_results.noise_realizations(successful_runs(i)).y_q_noise, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
        %     hold on;
        % end
        % xlabel('Sample');
        % ylabel('Angular Rate Noise [rad/s]');
        % title(sprintf('Angular Rate Noise Realizations (σ=%.4f)', noise.ang_rate_stand_dev));
        % grid on;
        
        % subplot(2, 2, 4);
        % for i = 1:run_example
        %     plot(mc_results.noise_realizations(successful_runs(i)).y_ax_noise, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
        %     hold on;
        % end
        % xlabel('Sample');
        % ylabel('Acceleration Noise [m/s²]');
        % title(sprintf('Acceleration Noise Realizations (σ≈%.4f)', noise.vel_stand_dev/sample_time));
        % grid on;
    end
    
    %% Save Results
    if save_results
        % Save full results
        filename = sprintf('monte_carlo_greybox_results_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
        save(filename, 'mc_results', 'mc_summary', 'noise', 'num_monte_carlo', 'real_parameters');
        
        if verbose
            fprintf('\nResults saved to: %s\n', filename);
        end
        mc_summary.saved_filename = filename;
    end
    
else
    if verbose
        fprintf('\nWarning: No successful identification runs. Consider:\n');
        fprintf('- Reducing noise levels in noise structure\n');
        fprintf('- Setting noise.Enabler = 0 for testing\n');
        fprintf('- Improving initial guess strategy\n');
        fprintf('- Checking model function implementation\n');
    end
    
    % Fill summary with empty/NaN values
    mc_summary.parameter_statistics = struct();
    mc_summary.error_statistics = struct();
end

end