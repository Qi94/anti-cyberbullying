clear;
clc;

load 'disucussion_experiment_1.mat'

% Create a figure with specified width and height
fig = figure('Position', [100, 100, 1200, 400]);

tiledlayout(1,2)

r_y_min = 0.1; % You've mentioned this value
r_y_max = 1.0;   % You've mentioned this value

x_scaled = linspace(r_y_min, r_y_max, experiment_times);

% First plot
nexttile
y_1 = profile_y_loss(1, :);
y_2 = profile_y_loss(2, :);
y_3 = profile_y_loss(3, :);

plot(x_scaled, y_1, '-bx', x_scaled, y_2, '-rx', x_scaled, y_3, '-gx', 'MarkerSize', 6)
l = legend('Facebook subnet', 'Youtube subnet', 'Twitter subnet');
%xlabel('Anti-Cyberbullying Turnover Rate');
xlabel('Anti-cyberbullying budget');
ylabel('{\it J_X(R^*_X,R^*_Y)}'); % <- y-axis label for the first plot
l.FontSize = 7;
%t.FontSize = 7;

% Second plot
nexttile
y_1 = profile_x_benifit(1, :);
y_2 = profile_x_benifit(2, :);
y_3 = profile_x_benifit(3, :);

plot(x_scaled, y_1, '-bx', x_scaled, y_2, '-rx', x_scaled, y_3, '-gx', 'MarkerSize', 6)
l = legend('Facebook subnet', 'Youtube subnet', 'Twitter subnet');
%xlabel('Anti-Cyberbullying Turnover Rate');
xlabel('Anti-cyberbullying budget');
ylabel('{\it J_Y(R^*_X,R^*_Y)}'); % <- y-axis label for the first plot
l.FontSize = 7;
%t.FontSize = 7;

% Save the figure as an .eps file
saveas(fig, 'asset/discusssion_experiment_1.eps', 'epsc');

clear;
clc;

load 'disucussion_experiment_2.mat'

% Create a figure with specified width and height
fig = figure('Position', [100, 100, 1200, 400]);

tiledlayout(1,2)

r_y_min = 0.001; % You've mentioned this value
r_y_max = 0.01;   % You've mentioned this value

x_scaled = linspace(r_y_min, r_y_max, experiment_times);

% First plot
nexttile
y_1 = profile_y_loss(1, :);
y_2 = profile_y_loss(2, :);
y_3 = profile_y_loss(3, :);

plot(x_scaled, y_1, '-bx', x_scaled, y_2, '-rx', x_scaled, y_3, '-gx', 'MarkerSize', 6)
l = legend('Facebook subnet', 'Youtube subnet', 'Twitter subnet');
%xlabel('Anti-Cyberbullying Turnover Rate');
xlabel('Anti-cyberbullying turnover rate');
ylabel('{\it J_X(R^*_X,R^*_Y)}'); % <- y-axis label for the first plot
l.FontSize = 7;
%t.FontSize = 7;

% Second plot
nexttile
y_1 = profile_x_benifit(1, :);
y_2 = profile_x_benifit(2, :);
y_3 = profile_x_benifit(3, :);

plot(x_scaled, y_1, '-bx', x_scaled, y_2, '-rx', x_scaled, y_3, '-gx', 'MarkerSize', 6)
l = legend('Facebook subnet', 'Youtube subnet', 'Twitter subnet');
%xlabel('Anti-Cyberbullying Turnover Rate');
xlabel('Anti-cyberbullying turnover rate');
ylabel('{\it J_Y(R^*_X,R^*_Y)}'); % <- y-axis label for the first plot
l.FontSize = 7;
%t.FontSize = 7;

% Save the figure as an .eps file
saveas(fig, 'asset/discusssion_experiment_2.eps', 'epsc');

clear;
clc;

load 'disucussion_experiment_3.mat'

% Create a figure with specified width and height
fig = figure('Position', [100, 100, 1200, 400]);

tiledlayout(1,2)

r_y_min = 0.002; % You've mentioned this value
r_y_max = 0.011;   % You've mentioned this value

x_scaled = linspace(r_y_min, r_y_max, experiment_times);

% First plot
nexttile
y_1 = profile_y_loss(1, :);
y_2 = profile_y_loss(2, :);
y_3 = profile_y_loss(3, :);

plot(x_scaled, y_1, '-bx', x_scaled, y_2, '-rx', x_scaled, y_3, '-gx', 'MarkerSize', 6)
l = legend('Facebook subnet', 'Youtube subnet', 'Twitter subnet');
%xlabel('Anti-Cyberbullying Turnover Rate');
xlabel('Cyberbullying turnover rate');
ylabel('{\it J_X(R^*_X,R^*_Y)}'); % <- y-axis label for the first plot
l.FontSize = 7;
%t.FontSize = 7;

% Second plot
nexttile
y_1 = profile_x_benifit(1, :);
y_2 = profile_x_benifit(2, :);
y_3 = profile_x_benifit(3, :);

plot(x_scaled, y_1, '-bx', x_scaled, y_2, '-rx', x_scaled, y_3, '-gx', 'MarkerSize', 6)
l = legend('Facebook subnet', 'Youtube subnet', 'Twitter subnet');
%xlabel('Anti-Cyberbullying Turnover Rate');
xlabel('Cyberbullying turnover rate');
ylabel('{\it J_Y(R^*_X,R^*_Y)}'); % <- y-axis label for the first plot
l.FontSize = 7;
%t.FontSize = 7;

% Save the figure as an .eps file
saveas(fig, 'asset/discusssion_experiment_3.eps', 'epsc');
