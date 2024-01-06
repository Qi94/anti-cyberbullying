% Load the necessary data
load plot-facebook10000.mat

% Create a figure with specified width and height
fig = figure('Position', [100, 100, 1200, 400]);

% Set up tiled layout
tiledlayout(1, 2)

% Plot xData
ax1 = nexttile; % Get handle to the axes
title('Cyberbullying Strategy');  
hold on
for i = 1:size(xData, 1)
    plot(time_points, xData(i, :), 'LineWidth', 1, 'DisplayName', sprintf('Loop %d', i));
end
% Highlight a region with a rectangle
rectangle('Position',[min(time_points) min(xData(:)) max(time_points)-min(time_points) max(xData(:))-min(xData(:))],'EdgeColor','r','LineStyle','--')
hold off
legend('Location', 'southwest');
xlabel('{\it t}');
ylabel('{\it R^*_x(t)}');
ax1.Box = 'on';

% Plot yData
ax2 = nexttile; % Get handle to the axes
title('Anti-cyberbullying Strategy');  
hold on
for i = 1:size(yData, 1)
    plot(time_points, yData(i, :), 'LineWidth', 1, 'DisplayName', sprintf('Loop %d', i));
end
% Highlight a region with a rectangle
rectangle('Position',[min(time_points) min(yData(:)) max(time_points)-min(time_points) max(yData(:))-min(yData(:))],'EdgeColor','r','LineStyle','--')
hold off
legend('Location', 'southwest');
xlabel('{\it t}');
ylabel('{\it R^*_y(t)}');
ax2.Box = 'on';

% Save the figure as an .eps file
saveas(fig, 'asset/Strategy_Facebook.eps', 'epsc');
