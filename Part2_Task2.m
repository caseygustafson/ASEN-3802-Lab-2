% Adam Cobb
% ASEN 3802 LAB 2 Part 2 Task 2
% 3/8/2026

clc; clear; close all;    % Clear Command Window, remove variables from workspace, and close all figure windows

% List of experimental data
files = [                 % Create a string array containing the names of the 5 experimental data files
    "Aluminum_25V_240mA"  % File 1: aluminum rod at 25 V, 240 mA
    "Aluminum_30V_290mA"  % File 2: aluminum rod at 30 V, 290 mA
    "Brass_25V_237mA"     % File 3: brass rod at 25 V, 237 mA
    "Brass_30V_285mA"     % File 4: brass rod at 30 V, 285 mA
    "Steel_22V_203mA"     % File 5: steel rod at 22 V, 203 mA
];

% Thermocouple positions measured from x0 [m]
x1 = 0.034925;           % Location of thermocouple 1 measured from x0, in meters
dx = 0.0127;             % Spacing between thermocouples, in meters
x = x1 + (0:7)*dx;       % Create vector of all 8 thermocouple positions: x1, x2, ..., x8

% Analytical steady-state slopes [°C/m]
H_an_vals = [91.086, 132.074, 101.679, 146.727, 544.059];
% Store the analytical steady-state slope H for each dataset, in the same order as 'files'

% Number of tail samples used for steady-state averaging
Nss = 125;               % Use the last 125 data points from each file to estimate steady-state temperatures

% Preallocate storage for results
T0_exp_list = zeros(length(files),1);   % Create column vector to store steady-state y-intercept T0 for each file
H_exp_list  = zeros(length(files),1);   % Create column vector to store experimental slope H_exp for each file


%% STEADY-STATE PROFILES
% This section computes steady-state averages and uses them to find T0 and H_exp

for i = 1:length(files)              % Loop through each of the 5 datasets

    % Reload dataset
    data = readmatrix(files(i));     % Read the current data file into a numeric matrix
    data(end-1:end,:) = [];          % Remove the last two rows (often extra junk/blank rows in exported lab data)
    data = data(all(~isnan(data),2),:);
    % Remove any row that contains NaN so only valid numeric data remains

    couples = data(:,2:9);           % Extract thermocouple data only (columns 2 through 9 = TC1 to TC8)

    % Steady-state tail region
    idx0 = max(1, size(couples,1) - Nss + 1);
    % Find the starting row index for the last Nss samples
    % max(1,...) prevents the index from going below 1 if the dataset is short

    % Compute steady-state averages and variation
    T_avg = mean(couples(idx0:end,:), 1);
    % Compute the mean temperature at each thermocouple over the steady-state tail region

    % Linear fit of temperatures
    p = polyfit(x, T_avg, 1);
    % Fit a straight line T(x) = p(1)*x + p(2) to the steady-state average temperature profile

    H_exp = p(1);                 % Experimental slope of the steady-state profile, in °C/m
    T0_exp = polyval(p, 0);       % Evaluate the fitted line at x = 0 to get the steady-state y-intercept T0

    % Analytical slope
    H_an = H_an_vals(i);          % Select the analytical steady-state slope for this file

    % Store results
    T0_exp_list(i) = T0_exp;      % Save the steady-state y-intercept from the experimental fit for this file
    H_exp_list(i)  = H_exp;       % Save the experimental slope for this file

end                              % End loop over files

%% RESULTS
% This section prints a summary table to the MATLAB Command Window

fprintf("%-20s  %-10s  %-12s  %-12s\n", "Dataset", "T0 [C]", "H_exp [C/m]", "H_an [C/m]");
% Print a formatted header row for the results table

for i = 1:length(files)          % Loop through each file again for printing
    fprintf("%-20s  %-10.3f  %-12.3f  %-12.3f\n", ...
        files(i), T0_exp_list(i), H_exp_list(i), H_an_vals(i));
    % Print dataset name, T0, H_exp, and H_an with 3 decimal places
end

%% Part 2 Task 2
% This section computes and plots the transient analytical solution (Model IA)
% and compares it to the experimental temperature data

% Thermal diffusivities [m^2/s]
alpha_vals = [ ...
    130/(2810*960), ...   % Aluminum thermal diffusivity alpha = k/(rho*cp)
    130/(2810*960), ...   % Same aluminum diffusivity for second aluminum case
    115/(8500*380), ...   % Brass thermal diffusivity
    115/(8500*380), ...   % Same brass diffusivity for second brass case
    16.2/(8000*500) ...   % Steel thermal diffusivity
];

% Full rod length [m]
L = 0.0127 * 7 + 0.0254 + 0.034925;
% Compute full rod length used in the analytical solution
% = 7 thermocouple spacings + extra 1 inch + distance from x0 to TC1

% Number of modes from Part 2 Task 1
N_modes = 6;             % Number of Fourier-series terms retained in the transient solution

for i = 1:length(files)  % Loop through all 5 datasets again

    % Load dataset  ← use files(i) directly, no .name needed
    data = readmatrix(files(i));   % Read current experimental data file
    data(end-1:end,:) = [];                    % Remove trailing junk rows
    data = data(all(~isnan(data),2),:);
    % Remove any rows containing NaN values

    % Experimental data
    t     = data(:,1);             % Extract time vector from column 1
    T_exp = data(:,2:9);           % Extract measured temperatures from TC1 through TC8

    % Model parameters
    H_an  = H_an_vals(i);          % Select analytical steady-state slope for current dataset
    T0    = T0_exp_list(i);        % Use steady-state y-intercept from Part 1 Task 1
    alpha = alpha_vals(i);         % Select thermal diffusivity for current material

    % Preallocate analytical temperatures
    T_model = zeros(length(t), 8);
    % Create matrix to store analytical temperatures for all times and all 8 thermocouples

    % Compute analytical solution at each thermocouple
    for j = 1:8                   % Loop over the 8 thermocouple locations
        xj = x(j);                % Store current thermocouple position

        % Steady-state term
        T_model(:,j) = T0 + H_an*xj;
        % Start the analytical solution with the steady-state contribution T_ss = T0 + H*x

        % Transient series
        for n = 1:N_modes         % Loop over the retained Fourier modes
            lambda_n = ((2*n - 1)*pi) / (2*L);
            % Compute eigenvalue lambda_n for the nth sine term

            b_n = (8 * H_an * L * (-1)^n) / (((2*n - 1)^2) * pi^2);
            % Compute Fourier coefficient b_n from your derived expression

            T_model(:,j) = T_model(:,j) + ...
                b_n * sin(lambda_n * xj) .* exp(-alpha * lambda_n^2 * t);
            % Add the nth transient term to the analytical temperature at this thermocouple
            % sin(lambda_n*xj) handles spatial dependence
            % exp(-alpha*lambda_n^2*t) handles decay with time
            % .* is element-by-element multiplication because t is a vector
        end
    end

    % ---------------- Plot ----------------
    figure;                       % Open a new figure window for the current dataset
    hold on;                      % Keep multiple lines on the same axes
    grid on;                      % Turn grid on for readability

    cmap = lines(8);              % Generate 8 distinct colors, one for each thermocouple
    h_an_leg  = gobjects(8,1);    % Preallocate graphics handles for analytical lines
    h_exp_leg = gobjects(8,1);    % Preallocate graphics handles for experimental lines

    for j = 1:8                   % Loop through all 8 thermocouples for plotting
        % Analytical = solid
        h_an_leg(j)  = plot(t, T_model(:,j), '-',  'Color', cmap(j,:), 'LineWidth', 1.8, ...
            'DisplayName', sprintf('TC%d Analytical  (x = %.1f mm)', j, x(j)*1000));
        % Plot analytical temperature for TC j as a solid line
        % Use color j from the colormap
        % Save the line handle in h_an_leg(j)
        % DisplayName creates the legend label

        % Experimental = dashed
        h_exp_leg(j) = plot(t, T_exp(:,j),   '--', 'Color', cmap(j,:), 'LineWidth', 1.0, ...
            'DisplayName', sprintf('TC%d Experimental (x = %.1f mm)', j, x(j)*1000));
        % Plot experimental temperature for TC j as a dashed line
        % Use the same color as the analytical line for that same thermocouple
        % Save the line handle in h_exp_leg(j)
    end

    xlabel('Time [s]');           % Label x-axis as time in seconds
    ylabel('Temperature [^\circ C]');
    % Label y-axis as temperature in degrees Celsius

    title(sprintf('Model IA Temperature vs Time — %s', strrep(files(i), '_', ' ')));
    % Create plot title
    % strrep replaces underscores in the file name with spaces for nicer formatting

    % Interleave analytical and experimental handles: TC1 An, TC1 Exp, TC2 An, TC2 Exp ...
    h_interleaved = reshape([h_an_leg, h_exp_leg]', [], 1);
    % Rearrange handles so legend order alternates:
    % TC1 Analytical, TC1 Experimental, TC2 Analytical, TC2 Experimental, etc.

    legend(h_interleaved, 'Location', 'eastoutside');
    % Create legend using the interleaved handles and place it outside the plot on the right

    text(0.02, 0.98, sprintf('H_{an} = %.3f ^\\circC/m\nT_0 = %.3f ^\\circC\nN = %d', ...
        H_an, T0, N_modes), ...
        'Units', 'normalized', 'VerticalAlignment', 'top');
    % Add text box inside the plot showing:
    % analytical slope H_an
    % steady-state intercept T0
    % number of modes used
    % normalized coordinates mean (0,0) is bottom-left and (1,1) is top-right of the axes

    hold off;                     % Release the current figure so the next loop starts fresh

end                               % End loop over all files
