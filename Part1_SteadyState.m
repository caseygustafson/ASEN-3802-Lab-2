clc; clear; close all;

% List of experimental data
files = [
    "Aluminum_25V_240mA"
    "Aluminum_30V_290mA"
    "Brass_25V_237mA"
    "Brass_30V_285mA"
    "Steel_22V_203mA"
];

% Thermocouple positions measured from x0 [m]
x1 = 0.034925;
dx = 0.0127;
x = x1 + (0:7)*dx;

% Analytical steady-state slopes[°C/m]
H_an_vals = [91.086, 132.074, 101.679, 146.727, 544.059];

% Number of tail samples used for steady-state averaging
Nss = 125;

% Preallocate storage for results
T0_exp_list = zeros(length(files),1);
H_exp_list  = zeros(length(files),1);

%% TEMPERATURE VS TIME
figure(1)

for i = 1:length(files)

    % Load dataset
    data = readmatrix(files(i));

    % RemoveNaN values
    data = data(all(~isnan(data),2),:);

    % Time and thermocouple temperatures
    t = data(:,1);
    couples = data(:,2:9);

    % Plot thermocouple temperatures vs time
    subplot(5,1,i)
    plot(t, couples)

    % Y-axis scaling with margin
    y_min = min(couples(:));
    y_max = max(couples(:));
    y_range = y_max - y_min;
    y_margin = max(1, 0.05*y_range);
    ylim([y_min - y_margin, y_max + y_margin])

    title(files(i) + " (Temperature of Thermocouples vs time)")
    xlabel("Time [s]")
    ylabel("Temp [C]")
    grid on
    
    % Legend
    if i == 1
        legend("TC1","TC2","TC3","TC4","TC5","TC6","TC7","TC8", ...
               "Location","eastoutside");
    end
end

%% STEADY-STATE PROFILES

figure(2)

for i = 1:length(files)

    % Reload dataset
    data = readmatrix(files(i));
    data(end-1:end,:) = [];
    data = data(all(~isnan(data),2),:);

    couples = data(:,2:9);

    % Steady-state tail region
    idx0 = max(1, size(couples,1) - Nss + 1);

    % Compute steady-state averages and variation
    T_avg = mean(couples(idx0:end,:), 1);
    T_std = std(couples(idx0:end,:), 0, 1);

    % Linear fit of temperatures
    p = polyfit(x, T_avg, 1);

    H_exp = p(1);                 % Experimental slope
    T0_exp = polyval(p, 0);       % Temperature at x0
    T_fit = polyval(p, x);        % Fitted profile

    % Analytical slope
    H_an = H_an_vals(i);

    % Starting analytical intercept at TC1
    T0_an = T_avg(1) - H_an*x(1);
    T_an  = T0_an + H_an*x;

    % Store results
    T0_exp_list(i) = T0_exp;
    H_exp_list(i)  = H_exp;

    % Plot steady-state data with fit and analytical model
    subplot(5,1,i)
    errorbar(x, T_avg, T_std, 'o')
    hold on
    plot(x, T_fit, '-')
    plot(x, T_an, '-')
    hold off

    % Y-axis scaling
    y_min = min([T_avg T_fit T_an]);
    y_max = max([T_avg T_fit T_an]);
    y_range = y_max - y_min;
    y_margin = max(1, 0.05*y_range);
    ylim([y_min - y_margin, y_max + y_margin])

    title(files(i) + " (Steady-state profile)")
    xlabel("x from x0 [m]")
    ylabel("Temp [C]")
    legend("Experimental","Polyfit","Analytical","Location","southeast")
    grid on
end

%% RESULTS

fprintf("\n%-20s  %-10s  %-12s  %-12s\n", "Dataset", "T0 [C]", "H_exp [C/m]", "H_an [C/m]");

for i = 1:length(files)
    fprintf("%-20s  %-10.3f  %-12.3f  %-12.3f\n", ...
        files(i), T0_exp_list(i), H_exp_list(i), H_an_vals(i));
end