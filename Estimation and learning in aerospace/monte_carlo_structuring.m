function [mc_results, mc_summary] = monte_carlo_structuring(noise, input_mc, ...
    t_3ord, Xu, Xq, Mu, Mq, Xd, Md, sample_time, fc, p_value, n_order, varargin)
%MONTE_CARLO_STRUCTURING Performs Monte Carlo simulation for PBSID parameter identification
%
% Syntax:
%   [mc_results, mc_summary] = monte_carlo_pbsid_identification(u_3ord, y_q, t_3ord, ...
%       Xu, Xq, Mu, Mq, Xd, Md, sample_time, fc, p_value, n_order)
%   [mc_results, mc_summary] = monte_carlo_pbsid_identification(..., 'Name', Value)
%
% Inputs:
%   u_3ord      - Control input vector
%   y_q         - Pitch rate measurements (single output for PBSID)
%   t_3ord      - Time vector
%   Xu,Xq,Mu,Mq,Xd,Md - True parameter values for reference
%   sample_time - Sampling time [s]
%   fc          - Cutting frequency for bettering solution
%   p_value     - Optimal p value for PBSID (from better_p function)
%   n_order     - System order detected by PBSID
%
% Optional Name-Value Pairs:
%   'num_monte_carlo' - Number of Monte Carlo runs (default: 100)
%   'noise_enabler'   - Enable/disable noise (default: 1)
%   'pos_noise_std'   - Position noise std dev [m] (default: 0.0011)
%   'vel_noise_std'   - Velocity noise std dev [m/s] (default: 0.01)
%   'att_noise_std'   - Attitude noise std dev [rad] (default: deg2rad(0.33))
%   'rate_noise_std'  - Angular rate noise std dev [rad/s] (default: deg2rad(1))
%   'model_function'  - Model function handle (default: @drone_model_PBSID)
%   'auto_p_selection' - Automatically determine p for each run (default: false)
%   'auto_order_selection' - Automatically determine order for each run (default: false)
%   'p_range'         - Range for automatic p selection [min_p, max_p] (default: [6, 40])
%   'verbose'         - Display progress (default: true)
%   'save_results'    - Save results to file (default: true)
%   'plot_results'    - Generate plots (default: true)
%
% Outputs:
%   mc_results  - Structure containing all Monte Carlo results
%   mc_summary  - Structure containing summary statistics
%
% Example:
%   [results, summary] = monte_carlo_structuring(u_input, pitch_rate, time_vec, ...
%       -0.5, 2.1, -1.2, -8.5, 0.2, -0.8, 0.01, 0.5, 12, 3, 'num_monte_carlo', 50);

%% Input Parsing
p = inputParser;
addRequired(p, 'Xu', @isnumeric);
addRequired(p, 'Xq', @isnumeric);
addRequired(p, 'Mu', @isnumeric);
addRequired(p, 'Mq', @isnumeric);
addRequired(p, 'Xd', @isnumeric);
addRequired(p, 'Md', @isnumeric);
addRequired(p, 'sample_time', @(x) isnumeric(x) && x > 0);
addRequired(p, 'fc', @(x) isnumeric(x) && x > 0);
addRequired(p, 'p_value', @(x) isnumeric(x) && x > 0);
addRequired(p, 'n_order', @(x) isnumeric(x) && x > 0);

addParameter(p, 'num_monte_carlo', 100, @(x) isnumeric(x) && x > 0);
addParameter(p, 'noise_enabler', 1, @(x) isnumeric(x) && (x == 0 || x == 1));
addParameter(p, 'pos_noise_std', 0.0011, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'vel_noise_std', 0.01, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'att_noise_std', deg2rad(0.33), @(x) isnumeric(x) && x >= 0);
addParameter(p, 'rate_noise_std', deg2rad(1), @(x) isnumeric(x) && x >= 0);
addParameter(p, 'model_function', @drone_model_PBSID, @(x) isa(x, 'function_handle'));
addParameter(p, 'auto_p_selection', false, @islogical);
addParameter(p, 'auto_order_selection', false, @islogical);
addParameter(p, 'p_range', [6, 40], @(x) isnumeric(x) && length(x) == 2);
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'save_results', true, @islogical);
addParameter(p, 'plot_results', true, @islogical);

parse(p, Xu, Xq, Mu, Mq, Xd, Md, sample_time, fc, p_value, n_order, varargin{:});

% Extract parsed values
num_monte_carlo = p.Results.num_monte_carlo;
model_fun = p.Results.model_function;
auto_p = p.Results.auto_p_selection;
auto_order = p.Results.auto_order_selection;
p_range = p.Results.p_range;
verbose = p.Results.verbose;
save_results = p.Results.save_results;
plot_results = p.Results.plot_results;

%% Noise Configuration
% default
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
mc_results.K_matrices = cell(num_monte_carlo, 1);  % Kalman gain for PBSID
mc_results.convergence_flags = zeros(num_monte_carlo, 1);
mc_results.p_values = zeros(num_monte_carlo, 1);  % Store p values used
mc_results.order_values = zeros(num_monte_carlo, 1);  % Store orders used
mc_results.noise_realizations = struct();
mc_results.pbsid_intermediate = cell(num_monte_carlo, 1);  % Store PBSID intermediate results

%% Display Configuration
if verbose
    fprintf('Starting Monte Carlo PBSID Parameter Identification...\n');
    fprintf('PBSID Configuration:\n');
    fprintf('  Cutting frequency: %.4f Hz\n', fc);
    fprintf('  P value: %d (auto: %s)\n', p_value, mat2str(auto_p));
    fprintf('  System order: %d (auto: %s)\n', n_order, mat2str(auto_order));
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
t=t_3ord;

A=[Xu, Xq, -9.81; Mu, Mq, 0; 0, 1, 0];

B=[Xd; Md; 0];

C=[1, 0, 0; 0, 1, 0; 0, 0, 1; Xu, Xq, 0]; 

D=[0; 0; 0; Xd];

simulation_time=t(end)-t(1);

assignin('base','ExcitationM', ExcitationM);
assignin('base','simulation_time', simulation_time);

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

    if mc_run == 1  % Only debug first run
    fprintf('Debug: ExcitationM size in workspace: [%d x %d]\n', size(ExcitationM));
    fprintf('Debug: simulation_time: %.3f\n', simulation_time);
    
    % Check if ExcitationM exists in base workspace
    if evalin('base', 'exist(''ExcitationM'', ''var'')')
        ExcitationM_check = evalin('base', 'ExcitationM');
        fprintf('ExcitationM in base workspace size: [%d x %d]\n', size(ExcitationM_check));
    else
        fprintf('ERROR: ExcitationM not found in base workspace!\n');
    end
end
    % --- Run simulation ---
    
    simulation_data = sim('Simulator_Single_Axis_MC','SrcWorkspace','current');

    % Extract & identify just like above
    u_3ord = simulation_data.Mtot; 
    y_q    = simulation_data.q;
    %y_ax   = simulation_data.ax;
    %% PBSID Identification Process with Noisy Data
    try
        % Input Assignment (PBSID uses single output)
        % y_3ord_noisy = y_q_noisy;
        Ts_3ord = sample_time;
        
        % Apply bettering solution with fixed fc
        %[u_3ord_filtered, y_3ord_filtered] = bettering_solution(u_3ord_noisy, y_3ord_noisy, fc, Ts_3ord);
        fs = 1/Ts_3ord; 
        [b,a] = butter(4, fc/(fs/2), 'low'); % Butterworth filter at 4 order
        u_filt = filtfilt(b, a, u_3ord);
        y_filt = filtfilt(b, a, y_q);

        % Normalization (porta i dati tra 0 e 1)
        u_norm = (u_filt - min(u_filt)) / (max(u_filt) - min(u_filt));
        y_norm = (y_filt - min(y_filt)) / (max(y_filt) - min(y_filt));

        % Frame length of Savitzky-Golay filter
        framelen = fs/fc * 5; % 
        framelen = round(framelen) + (mod(round(framelen), 2) == 0) * (2 * (round(framelen) < framelen) - 1);
        framelen = 21;
        % Smoothing with Savitzky-Golay filter
        u_3ord_filtered = sgolayfilt(u_norm, 3, framelen);
        y_3ord_filtered = sgolayfilt(y_norm, 3, framelen);
        % Determine p value for current run
        % if auto_p
        %     % Automatically determine optimal p for this noisy realization
        %     min_p = p_range(1); 
        %     max_p = p_range(2);
        %     current_p = better_p(u_3ord_filtered, y_3ord_filtered, min_p, max_p, 3, Ts_3ord, t_3ord);
        % else
            current_p = p_value;
        % end
        mc_results.p_values(mc_run) = current_p;
        
        % PBSID Part 1
        [D_PBSID, S_PBSID, V_PBSID, Y, N] = pbsid_1part(u_3ord_filtered, y_3ord_filtered(:,1), current_p, 1);
        
        % Determine system order for current run
        if auto_order
            % Automatically determine order (could implement SVD-based selection)
            % For now, use a simple heuristic or the provided order with some variation
            current_n = max(1, n_order + randi([-1, 1]));  % Small random variation
        else
            current_n = n_order;
        end
        mc_results.order_values(mc_run) = current_n;
        
        % PBSID Part 2
        [A_PBSID, B_PBSID, C_PBSID, K_PBSID] = pbsid_2part(D_PBSID, current_n, S_PBSID, V_PBSID, Y, N, current_p, u_3ord_filtered, y_3ord_filtered);
        
        % Set D matrix (usually zero for PBSID)
        if size(D_PBSID, 1) == 0 || size(D_PBSID, 2) == 0
            D_PBSID = zeros(size(C_PBSID, 1), size(B_PBSID, 2));
        end
        
        % Generate intermediate PBSID output for greyest structuring
        y_id = lsim(ss(A_PBSID, B_PBSID, C_PBSID, D_PBSID), u_3ord_filtered, t_3ord);
        
        % Initial parameter guess with random perturbation
        theta0 = [Xu Xq Mu Mq Xd Md] + 0.01 * (-0.5 + rand(1,6)).*[Xu Xq Mu Mq Xd Md];
        
        % Greyest structuring
        [identification_PBSID, error_PBSID] = greyest_structuring(u_3ord_filtered, y_id, Ts_3ord, theta0, model_fun, real_parameters);
        
        % Extract final structured matrices
        [A_final, B_final, C_final, D_final] = drone_model_grey(identification_PBSID.parameters, 0);
        
        % Store results
        mc_results.identified_params(mc_run, :) = identification_PBSID.parameters';
        mc_results.errors(mc_run, :) = error_PBSID;
        mc_results.A_matrices{mc_run} = A_final;
        mc_results.B_matrices{mc_run} = B_final;
        mc_results.C_matrices{mc_run} = C_final;
        mc_results.D_matrices{mc_run} = D_final;
        mc_results.K_matrices{mc_run} = K_PBSID;
        mc_results.convergence_flags(mc_run) = 1;  % Successful identification
        
        % Store intermediate PBSID results
        mc_results.pbsid_intermediate{mc_run} = struct(...
            'A_PBSID', A_PBSID, 'B_PBSID', B_PBSID, ...
            'C_PBSID', C_PBSID, 'D_PBSID', D_PBSID, ...
            'y_id', y_id);
        
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
        mc_results.K_matrices{mc_run} = [];
        mc_results.convergence_flags(mc_run) = 0;  % Failed identification
        mc_results.p_values(mc_run) = NaN;
        mc_results.order_values(mc_run) = NaN;
        mc_results.pbsid_intermediate{mc_run} = [];
    end
    
    %% Store noise realizations for this run
    % mc_results.noise_realizations(mc_run).y_q_noise = ang_rate_noise_q;
    % mc_results.noise_realizations(mc_run).input_noise = input_noise;
end

if verbose
    elapsed_time = toc;
    fprintf('\nMonte Carlo PBSID simulation completed in %.2f seconds\n', elapsed_time);
else
    elapsed_time = 0;
end

%% Statistical Analysis
successful_runs = find(mc_results.convergence_flags == 1);
num_successful = length(successful_runs);
success_rate = num_successful / num_monte_carlo * 100;

if verbose
    fprintf('\nMonte Carlo PBSID Results Summary:\n');
    fprintf('Success rate: %.1f%% (%d/%d runs)\n', success_rate, num_successful, num_monte_carlo);
    
    if auto_p && num_successful > 0
        successful_p = mc_results.p_values(successful_runs);
        fprintf('P value statistics: mean=%.1f, std=%.2f, range=[%d, %d]\n', ...
            mean(successful_p), std(successful_p), min(successful_p), max(successful_p));
    end
    
    if auto_order && num_successful > 0
        successful_orders = mc_results.order_values(successful_runs);
        fprintf('Order statistics: mean=%.1f, std=%.2f, range=[%d, %d]\n', ...
            mean(successful_orders), std(successful_orders), min(successful_orders), max(successful_orders));
    end
end

%% Initialize summary structure
mc_summary = struct();
mc_summary.num_runs = num_monte_carlo;
mc_summary.noise_config = noise;
mc_summary.pbsid_config = struct('fc', fc, 'p_value', p_value, 'n_order', n_order, ...
                                 'auto_p', auto_p, 'auto_order', auto_order);
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
        figure('Name', 'Monte Carlo PBSID Parameter Distributions', 'Position', [100 100 1200 800]);
        
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
        
        % PBSID-specific plots
        figure('Name', 'Monte Carlo PBSID Analysis', 'Position', [200 200 1000 800]);
        
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
        
        % P value and order distribution (if auto selection enabled)
        if auto_p
            subplot(2, 2, 3);
            histogram(mc_results.p_values(successful_runs), 15);
            xlabel('P Value');
            ylabel('Frequency');
            title('P Value Distribution (Auto Selection)');
            grid on;
        % else
        %     subplot(2, 2, 3);
        %     run_example = min(10, num_successful);
        %     for i = 1:run_example
        %         plot(mc_results.noise_realizations(successful_runs(i)).y_q_noise, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
        %         hold on;
        %     end
        %     xlabel('Sample');
        %     ylabel('Angular Rate Noise [rad/s]');
        %     title(sprintf('Noise Realizations (σ=%.4f)', noise.ang_rate_stand_dev));
        %     grid on;
        end
        
        if auto_order
            subplot(2, 2, 4);
            histogram(mc_results.order_values(successful_runs), 10);
            xlabel('System Order');
            ylabel('Frequency');
            title('Order Distribution (Auto Selection)');
            grid on;
        else
            subplot(2, 2, 4);
            % Show example of PBSID intermediate results
            if ~isempty(mc_results.pbsid_intermediate{successful_runs(1)})
                example_y_id = mc_results.pbsid_intermediate{successful_runs(1)}.y_id;
                plot(t_3ord, example_y_id, 'b-', 'LineWidth', 1);
                xlabel('Time [s]');
                ylabel('PBSID Intermediate Output');
                title('Example PBSID Intermediate Result');
                grid on;
            end
        end
    end
    
    %% Save Results
    if save_results
        filename = sprintf('monte_carlo_pbsid_results_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
        save(filename, 'mc_results', 'mc_summary', 'noise', 'num_monte_carlo', 'real_parameters');
        
        if verbose
            fprintf('\nResults saved to: %s\n', filename);
        end
        mc_summary.saved_filename = filename;
    end
    
else
    if verbose
        fprintf('\nWarning: No successful PBSID identification runs. Consider:\n');
        fprintf('- Reducing noise levels\n');
        fprintf('- Adjusting cutting frequency (fc)\n');
        fprintf('- Modifying p_value or enabling auto_p_selection\n');
        fprintf('- Checking system order selection\n');
        fprintf('- Verifying PBSID function implementations\n');
    end
    
    % Fill summary with empty/NaN values
    mc_summary.parameter_statistics = struct();
    mc_summary.error_statistics = struct();
end

end

