%{
Author: QC
CreateOn: 24/09/2022 

Reproduce of paper 'Game-theoretic modelling and analysis of cyberbullying spreading on OSNs' 
%}

clear;
clc;

%% ===============parameter initiation===============

% ---iteration parameters---
% number of iterations
n = 10^3;
% convergence error
epsilon = 10^-6;

% Initialize the execution times array
iterationDurations = zeros(n, 1);

% Record the start time of the entire execution
totalStartTime = tic;


% ---model parameters---
% control period T
T = 5;

% ---load network structure---
% load the adjacency matrix
%load -ASCII asset/twitter100.mat
load('asset/bingHam10000.mat');

% adjacency matrix a
%a = twitter100;
a = matrix;

% number of individuals in the net
N = size(a,1);

% upper bound on instant increase rate of conducting cyberbullying strategy at time
r_x_max = 1;

% upper bound on instant increase rate of conducting anti-cyberbullying strategy at time
r_y_max = 1;

% likelihood when a by-standing user (anti-cyberbullying user, respectively) becomes cyberbullying because of the influence of a cyberbullying user connected with him/her on OSN
beta = [0.012, 0.014];

% likelihood when a by-standing user (cyberbullying user, respectively) becomes anti-cyberbullying because of the influence of a anti-cyberbullying user connected with him/her on OSN
gamma = [0.013, 0.015];

% likelihood when a by-standing user becomes cyberbullying because of the influence of a anti-cyberbullying user connected with him/her on OSN
zeta = 0.0023;

% likelihood when a by-standing user becomes anti-cyberbullying because of the influence of a cyberbullying user connected with him/her on OSN
eta = 0.0025;

% likelihood when a by-standing user (anti-cyberbullying user, respectively) becomes cyberbullying because of the influence of R_X
% type of function f. 'convex', 'concave' or 'normal'
f_a = @SCR_f_a_normal;
f_b = @SCR_f_b_normal;

% likelihood when a by-standing user (cyberbullying user, respectively) becomes anti-cyberbullying because of the influence of R_Y
% type of function g. 'convex', 'concave' or 'normal'
g_c = @SCR_g_c_normal;
g_b = @SCR_g_b_normal;

% likelihood when a by-standing user becomes anti-cyberbullying because of the influence of R_X
v = @SCR_v_normal;

% likelihood when a by-standing user becomes cyberbullying because of the influence of R_Y
m = @SCR_m_normal;

% average benefit per unit time of X owing to a cyberbullying user (w_X), and the average loss per unit time of $Y$ owing to a cyberbullying user (w_Y)
w = [0.00232, 0.00233];

% The initial state of the OSN at the beginning, randomly generated N * M martrix
rng(1);
O = rand(N, 3);
rowsum = sum(O,2);
O = bsxfun(@rdivide, O, rowsum);

%% ===============variable initialization===============
% differential time
dt = T/n;
% discretize the time period
time_points = 0:dt:T;
% network state E
E_C = zeros(N, n+1);
E_A = zeros(N, n+1);
% costate lambdaA
lambdaAi = zeros(N, n+1);
% costate lambdc
lambdaCi = zeros(N, n+1);
% costate mua
muAi = zeros(N, n+1);
% costate muc
muCi = zeros(N, n+1);
% initial control <- all zero control
r_x = zeros(1, n+1);
r_y = ones(1, n+1);

% Initialize xData and yData with zero rows, as m is unknown initially
xData = zeros(0, n+1);
yData = zeros(0, n+1);

% Initialize loopCounter
loopCounter = 1;

%% ===============forward-backward sweep===============
% ---initial state---
% E(0) = E_0
E_C(:,1) = O(:,1);
E_A(:,1) = O(:,2);

% Initialize mainLoop counter
mainLoop = 1;

% ---forward-backward Euler algorithm---
% number of iteration
iter = 0;
while true
    % Start timing this iteration
    iterationStartTime = tic;

    % --preliminary preparation--
    %report the number of iterations
    iter = iter+1;
    disp(['progressing ',num2str(iter),'th iteration...']);
    % retain u of the previous iteration
    % to determine if the algorithm terminated
    r_x_prev = r_x;
    r_y_prev = r_y;

    % --forward calculate state E--
    disp('  forward sweep');
    % A(t+1) = A(t) + dt * (dAi+dCi)
    for t = 1:n
        % calculate each node separately for simplicity
        dc = (SCR_f_b_normal(r_x(t)) + beta(2) * a * E_C(:,t) + SCR_m_normal(r_y(t)) + zeta * a * E_A(:,t)) .* ...
            (1 - E_C(:,t) - E_A(:,t)) + (SCR_f_a_normal(r_x(t)) + beta(1) * a * E_C(:,t)) .* E_A(:,t) - ...
            (SCR_g_c_normal(r_y(t)) + gamma(2) * a * E_A(:,t)) .* E_C(:,t);
        
        da = (SCR_g_b_normal(r_y(t)) + gamma(1) * a * E_A(:,t) + SCR_v_normal(r_x(t)) + eta * a * E_C(:,t)) .* ...
            (1 - E_C(:,t) - E_A(:,t)) + (SCR_g_c_normal(r_y(t)) + gamma(2) * a * E_A(:,t)) .* E_C(:,t) - ...
            (SCR_f_a_normal(r_x(t)) + beta(1) * a * E_C(:,t)) .* E_A(:,t);
        
        % update
        E_C(:,t+1) = E_C(:,t) + dt * dc;
        E_A(:,t+1) = E_A(:,t) + dt * da;
        
        % make sure E_A and E_C is in [0,1]
        E_C(:,t+1) = state_adjust(E_C(:,t+1), 'c');
        E_A(:,t+1) = state_adjust(E_A(:,t+1), 'a');  

    end % end forward sweep
    
    
    % --backward calculate costate lambda--
    disp('  backward sweep');
    % lambda(t-1) = lambda(t) - dt * (dlambdaAi+dlambdaCi)
    % mu(t-1) = mu(t) - dt * (dmuAi+dmuCi)
    for t = n+1:-1:2
        % calculate the result
        % TODO: w_x = w(1), w_y = w(2) ?
        dlambda_C = -w(1) + (SCR_f_b_normal(r_x(t)) + SCR_g_c_normal(r_y(t)) + SCR_m_normal(r_y(t)) + ...
            beta(2) * a * E_C(:,t) + zeta * a * E_A(:,t) + gamma(2) * a * E_A(:,t)) .* lambdaCi(:,t) - ...
            a * ((beta(2) * lambdaCi(:,t) + eta * lambdaAi(:,t)) .* (1 - E_A(:,t) - E_C(:,t)) - ...
            beta(1) * E_A(:,t) .* (lambdaAi(:,t) - lambdaCi(:,t))) + ...
            (SCR_g_b_normal(r_y(t)) + SCR_v_normal(r_x(t)) - SCR_g_c_normal(r_y(t)) + ...
            eta * a * E_C(:,t) + gamma(1) * a * E_A(:,t) - gamma(2) * a * E_A(:,t)) .* lambdaAi(:,t);
        
        dlambda_A = (SCR_f_b_normal(r_x(t)) + SCR_m_normal(r_y(t)) - SCR_f_a_normal(r_x(t)) + ...
            (beta(2) - beta(1)) * a * E_C(:,t) + zeta * a * E_A(:,t)) .* lambdaCi(:,t) + ...
            gamma(2) * a * E_C(:,t) .* (lambdaCi(:,t) - lambdaAi(:,t)) + ...
            (SCR_f_a_normal(r_x(t)) + SCR_g_b_normal(r_y(t)) + SCR_v_normal(r_x(t)) + ...
            (beta(1) + eta) * a * E_C(:,t) + gamma(1) * a *E_A(:,t)) .* lambdaAi(:,t) - ...
            a * ((gamma(1) * lambdaAi(:,t) + zeta * lambdaCi(:,t)) .* (1 - E_A(:,t) - E_C(:,t)));
        
        dmu_C = -w(2) + (SCR_f_b_normal(r_x(t)) + SCR_g_c_normal(r_y(t)) + SCR_m_normal(r_y(t)) + ...
            beta(2) * a * E_C(:,t) + zeta * a * E_A(:,t) + gamma(2) * a * E_A(:,t)) .* muCi(:,t) - ...
            a * ((beta(2) * muCi(:,t) + eta* muAi(:,t)) .* (1 - E_A(:,t) - E_C(:,t)) - ...
            beta(1) * E_A(:,t) .* (muAi(:,t) - muCi(:,t))) + ...
            (SCR_g_b_normal(r_y(t)) + SCR_v_normal(r_x(t)) - SCR_g_c_normal(r_y(t)) + ...
            eta * a * E_C(:,t) + gamma(1) * a * E_A(:,t) - gamma(2) * a * E_A(:,t)) .* muAi(:,t);
        
        dmu_A = (SCR_f_b_normal(r_x(t)) + SCR_m_normal(r_y(t)) - SCR_f_a_normal(r_x(t)) + ...
            (beta(2) - beta(1)) * a * E_C(:,t) + zeta * a * E_A(:,t)) .* muCi(:,t) + ...
            gamma(2) * a * E_C(:,t) .* (muCi(:,t) - muAi(:,t)) + ...
            (SCR_f_a_normal(r_x(t)) + SCR_g_b_normal(r_y(t)) + SCR_v_normal(r_x(t)) + ...
            (beta(1) + eta) * a * E_C(:,t) + gamma(1) * a * E_A(:,t)) .* muAi(:,t) - ...
            a * ((gamma(1) * muAi(:,t) + zeta * muCi(:,t)) .* (1 - E_A(:,t) - E_C(:,t)));
        
        % update
        lambdaCi(:,t-1) = lambdaCi(:,t) - dt * dlambda_C;
        lambdaAi(:,t-1) = lambdaAi(:,t) - dt * dlambda_A;
        muCi(:,t-1) = muCi(:,t) - dt * dmu_C;
        muAi(:,t-1) = muAi(:,t) - dt * dmu_A;
    end
    
    % --update r_x and r_y--
    disp('  update control');
    % using different methods in different situations
    % using fminbnd() to find the maximum

    for t = 1:n+1
        term_final_x = @(x) -1 * ( f_b(x) * ((1-E_C(:,t)-E_A(:,t))' * lambdaCi(:,t)) + ...
            v(x) * ((1 - E_C(:,t) - E_A(:,t))' * lambdaAi(:,t)) + ...
            f_a(x) * (E_A(:,t)' * (lambdaCi(:,t) - lambdaAi(:,t))) - x );
        r_x(:,t) = fminbnd(term_final_x,0,r_x_max);
        
    end

    % Store r_x at the end of each iteration
    xData(loopCounter, :) = r_x;

    % using fminbnd() to find the minimum
    for t = 1:n+1
        term_final_y = @(y) m(y) * ((1 - E_C(:,t) - E_A(:,t))' * muCi(:,t)) + ...
            g_b(y) * ((1 - E_C(:,t) - E_A(:,t))' * muAi(:,t)) + ...
            g_c(y) * (E_C(:,t)' * (muAi(:,t) - muCi(:,t))) + y;
        r_y(:,t) = fminbnd(term_final_y,0,r_y_max);

    end

    % Store r_y at the end of each iteration
    yData(loopCounter, :) = r_y;

    % Increment loopCounter for the next iteration
    loopCounter = loopCounter + 1;

    % Store the duration of this iteration
    iterationDurations(iter) = toc(iterationStartTime);

    % --determine if the algorithm reaches a convergence--
    % if the difference between r_x_prev (r_y_prev) and r_x (r_y) is less than epsilon then the algorithm reaches a convergence
    dif_r_x = abs(r_x_prev - r_x);
    max_dif_r_x = max(dif_r_x,[],'all');
    disp(['  sum of difference between R(X) 2 iteration: ',num2str(sum(dif_r_x,'all'))]);
    disp(['  max difference between R(X) 2 iteration: ',num2str(max_dif_r_x)]);

    dif_r_y = abs(r_y_prev - r_y);
    max_dif_r_y = max(dif_r_y,[],'all');
    disp(['  sum of difference between R(Y) 2 iteration: ',num2str(sum(dif_r_y,'all'))]);
    disp(['  max difference between R(Y) 2 iteration: ',num2str(max_dif_r_y)]);

    if (max_dif_r_x < epsilon) & (max_dif_r_y < epsilon)
        disp(['--algorithm terminated at ',num2str(iter),'th iteration--']);
        break;
    end
end

% Record the total execution time
totalExecutionTime = toc(totalStartTime);

% Save the iteration durations and total execution time to the workspace
%save('execution_times_facebook1k.mat', 'iterationDurations', 'totalExecutionTime');


filename = 'plot-twitter10000.mat';
%filename = 'plot-Facebook1000.mat';
save(filename)

%% ===============result analysis===============
% ---plot control r_x---
%hold on
%figure;
%subplot(time_points,r_x(1,:),'r');
% ---plot control r_y---
%figure;
%subplot(time_points,r_y(1,:),'g');

% tiledlayout(1,2)
% nexttile
% plot(time_points,r_x(1,:),'r');
% 
% nexttile
% plot(time_points,r_y(1,:),'g');

