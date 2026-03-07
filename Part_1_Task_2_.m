% Adam Cobb
% ASEN 3802 LAB 2 Part 1 Task 2
% 2/25/2026


clear; 
clc; 
close all;

files = dir('*mA');   % Creates a file directory
inch_to_meters  = 0.0254;       % m/in

dx_in = 0.5; % change of spacing of thermocouples in inches
dx    = dx_in * inch_to_meters; % inches to meters for change in spacing

x1 = (11/8) * inch_to_meters; % Th1 location from x0, in meters
x  = x1 + (0:7) * dx; % 8 thermocouples

% loops through and runs and computes slopes
for i = 1:length(files)
    fname = files(i).name;

    % Read file (skip header row)
    data = readmatrix(fname, 'NumHeaderLines', 1);

    % Columns: time, CH1..CH8
    % Initial temps (at first row; typically t=0)
    Tinit = data(1, 2:9)';   % 8x1 vector (Th1..Th8)

    % Linear fit: T = Mexp*x + b
    p = polyfit(x, Tinit, 1); % p(1)=slope, p(2)=intercept
    Mexp = p(1);  % grab the slope in °C/m

    % Predicted temps at thermocouple locations
    Tfit = polyval(p, x);

    % plot overlay fit onto raw data
    figure; hold on; grid on;
    plot(x, Tinit, 'o', 'LineWidth', 1.5, 'DisplayName', 'Raw initial data');
    plot(x, Tfit,  '-', 'LineWidth', 2.0, 'DisplayName', 'Linear fit (predicted)');

    % show the predicted values as markers too
    plot(x, Tfit, 'kx', 'LineWidth', 1.5, 'DisplayName', 'Predicted at TCs');

    xlabel('x (m)');
    ylabel('T (^{\circ}C)');
    
    [~, base, ~] = fileparts(fname);
    title(sprintf('Initial Temperature Distribution — %s', strrep(base, '_', ' ')));
    
    text(0.05, 0.95, sprintf('M_{exp} = %.4f °C/m', Mexp), ...
    'Units', 'normalized', 'FontSize', 11, 'VerticalAlignment', 'top');

    legend('Location','southeast');

    % Print slope in command window
    fprintf('%s  -->  Mexp = %.6f (°C/m)\n', fname, Mexp);

    % Save the figure as a PNG file
    saveas(gcf, sprintf('%s_initial_temp_distribution.png', strrep(base, ' ', '_')));

end
