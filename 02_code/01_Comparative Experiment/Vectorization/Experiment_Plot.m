% Load the necessary data
clear;
clc;


load 'Experiment-cyberbullying-Youtube100.mat';

% Create a figure with specified width and height
fig = figure('Position', [100, 100, 1200, 400]);

tiledlayout(1,2)

% First plot for J_X Benefit
nexttile
x_c = 1:experiment_times+1;
y_c = [J_y_optimality, J_X];
stem(x_c,y_c, 'filled')
hold on
stem(x_c([1]), y_c([1]), 'r','MarkerFaceColor','r');
hold off
ylim([0 65])
xlim([0 100])
title('J_X Benefit')
xlabel('Anti-cyberbullying strategies');
ylabel('{\it J_X(R_X,R^*_Y)}'); % <- y-axis label for the first plot

load 'Experiment-anti-cyberbullying-Youtube100.mat';

% Second plot for J_Y Loss
nexttile
x_a = 1:experiment_times+1;
y_a = [J_x_optimality, J_Y];
stem(x_a,y_a, 'filled');
hold on
stem(x_a([1]), y_a([1]), 'r','MarkerFaceColor','r');
hold off
ylim([0 65])
xlim([0 100])
title('J_Y Loss')
xlabel('Cyberbullying strategies');
ylabel('{\it J_Y(R^*_X,R_Y)}'); % <- y-axis label for the second plot

% Save the figure as an .eps file
saveas(fig, 'asset/Evaluation-YoutubeSubnet.eps', 'epsc');


% Load the necessary data
clear;
clc;


load 'Experiment-cyberbullying-facebook100.mat';

% Create a figure with specified width and height
fig = figure('Position', [100, 100, 1200, 400]);

tiledlayout(1,2)

% First plot for J_X Benefit
nexttile
x_c = 1:experiment_times+1;
y_c = [J_y_optimality, J_X];
stem(x_c,y_c, 'filled')
hold on
stem(x_c([1]), y_c([1]), 'r','MarkerFaceColor','r');
hold off
ylim([0 65])
xlim([0 100])
title('J_X Benefit')
xlabel('Anti-cyberbullying strategies');
ylabel('{\it J_X(R_X,R^*_Y)}'); % <- y-axis label for the first plot

load 'Experiment-anti-cyberbullying-facebook100.mat';

% Second plot for J_Y Loss
nexttile
x_a = 1:experiment_times+1;
y_a = [J_x_optimality, J_Y];
stem(x_a,y_a, 'filled');
hold on
stem(x_a([1]), y_a([1]), 'r','MarkerFaceColor','r');
hold off
ylim([0 65])
xlim([0 100])
title('J_Y Loss')
xlabel('Cyberbullying strategies');
ylabel('{\it J_Y(R^*_X,R_Y)}'); % <- y-axis label for the second plot

% Save the figure as an .eps file
saveas(fig, 'asset/Evaluation-FacebookSubnet.eps', 'epsc');

load 'Experiment-cyberbullying-Twitter100.mat';

% Create a figure with specified width and height
fig = figure('Position', [100, 100, 1200, 400]);

tiledlayout(1,2)

% First plot for J_X Benefit
nexttile
x_c = 1:experiment_times+1;
y_c = [J_y_optimality, J_X];
stem(x_c,y_c, 'filled')
hold on
stem(x_c([1]), y_c([1]), 'r','MarkerFaceColor','r');
hold off
ylim([0 65])
xlim([0 100])
title('J_X Benefit')
xlabel('Anti-cyberbullying strategies');
ylabel('{\it J_X(R_X,R^*_Y)}'); % <- y-axis label for the first plot

load 'Experiment-anti-cyberbullying-Twitter100.mat';

% Second plot for J_Y Loss
nexttile
x_a = 1:experiment_times+1;
y_a = [J_x_optimality, J_Y];
stem(x_a,y_a, 'filled');
hold on
stem(x_a([1]), y_a([1]), 'r','MarkerFaceColor','r');
hold off
ylim([0 65])
xlim([0 100])
title('J_Y Loss')
xlabel('Cyberbullying strategies');
ylabel('{\it J_Y(R^*_X,R_Y)}'); % <- y-axis label for the second plot

% Save the figure as an .eps file
saveas(fig, 'asset/Evaluation-TwitterSubnet.eps', 'epsc');
