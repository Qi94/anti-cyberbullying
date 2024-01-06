%tiledlayout(1,3)
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Load the Facebook network data
load('asset/cal_10000.mat');
A_facebook100 = matrix;
G_facebook100 = graph(A_facebook100,'lower');
nexttile
p_facebook100 = plot(G_facebook100, 'Layout', 'force', 'EdgeAlpha', 0.1, 'NodeColor', 'r', 'NodeLabel', []);
p_facebook100.LineStyle = "-";
p_facebook100.Marker = ".";
p_facebook100.MarkerSize = 2;  % Smaller marker size
[t,s] = title('10000-node Facebook Network');
t.FontSize = 15;

% Load the Twitter network data
load('asset/bingHam10000.mat');
A_twitter100 = matrix;
G_twitter100 = graph(A_twitter100,'lower');
nexttile
p_twitter100 = plot(G_twitter100, 'Layout', 'force', 'EdgeAlpha', 0.1, 'NodeColor', 'r', 'NodeLabel', []);
p_twitter100.LineStyle = "-";
p_twitter100.Marker = ".";
p_twitter100.MarkerSize = 2;  % Smaller marker size
[t,s] = title('10000-node Twitter Network');
t.FontSize = 15;

% Load the Youtube network data
load('asset/duke_10000.mat');
A_youtube100 = matrix;
G_youtube100 = graph(A_youtube100,'lower');
nexttile
p_youtube100 = plot(G_youtube100, 'Layout', 'force', 'EdgeAlpha', 0.1, 'NodeColor', 'r', 'NodeLabel', []);
p_youtube100.LineStyle = "-";
p_youtube100.Marker = ".";
p_youtube100.MarkerSize = 2;  % Smaller marker size
[t,s] =title('10000-node Youtube Network');
t.FontSize = 15;
