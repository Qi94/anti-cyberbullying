
clear;
clc;

% ---load network structure---
Number_of_Network = 3;

%load -ASCII asset/Facebook100.mat
load('asset/cal_10000.mat');
Network_F = matrix;

%load -ASCII asset/Youtube100.mat
load('asset/duke_10000.mat');
Network_Y = matrix;

%load -ASCII asset/Twitter100.mat
load('asset/bingHam10000.mat');
Network_T = matrix;




% upper bound on instant increase rate of conducting anti-cyberbullying strategy at time
experiment_times = 10;
r_y_min = 0.1;
frequency = 0.1;

r_y_max_list = zeros(1, experiment_times);

profile_y_loss = zeros(Number_of_Network, experiment_times);
profile_x_benifit = zeros(Number_of_Network, experiment_times);

for n = 1:experiment_times
    r_y_max_list(1, n) = r_y_min;
    r_y_min = r_y_min + frequency;
end


% for n = 1:experiment_times
%     x = r_y_max_list(1, n);
%     profile_y_loss(1, n) = x;
% end

%{
Author: QC
CreateOn: 24/09/2022 

Reproduce of paper 'Game-theoretic modelling and analysis of cyberbullying spreading on OSNs' 
%}

%% ===============parameter initiation===============
% ---iteration parameters---
% number of iterations
n = 10^3;
% convergence error
epsilon = 10^-6;

% ---model parameters---
% control period T
T = 5;

% upper bound on instant increase rate of conducting cyberbullying strategy at time
r_x_max = 1;


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

for network_iter = 1:Number_of_Network
    if network_iter == 1
        a = Network_F;
    elseif network_iter == 2 
        a = Network_Y;
    elseif network_iter == 3
        a = Network_T;
    end
    % number of individuals in the net
    N = size(a,1);
    for experiment_iter = 1:experiment_times
        % The initial state of the OSN at the beginning, randomly generated N * M martrix
        rng(1);
        O = rand(N, 3);
        rowsum = sum(O,2);
        O = bsxfun(@rdivide, O, rowsum);
        %% ===============variable initialization===============
        % upper bound on instant increase rate of conducting anti-cyberbullying strategy at time
        r_y_max = r_y_max_list(1, experiment_iter);
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
                
        
        %% ===============Total Benefit & Loss Measurement===============
        % measure the total cyberbullying benefit and anti-cyberbullyong loss
        % when using optimality startegy profiles
        J_x_optimality = 0;
        J_y_optimality = 0;
        
        for t = 1:n
            term_J_X = 0;
            term_J_Y = 0;    
            for i = 1:N
                term_J_X = term_J_X + w(1)*E_C(i,t);      
                term_J_Y = term_J_Y + w(2)*E_C(i,t);  
            end
            J_x_optimality = J_x_optimality + term_J_X*dt - r_x(1,t)*dt;
            J_y_optimality = J_y_optimality + term_J_Y*dt + r_y(1,t)*dt;
        end
        
        profile_x_benifit(network_iter, experiment_iter) = J_x_optimality;
        profile_y_loss(network_iter, experiment_iter) = J_y_optimality;
    
    end
end

filename = 'disucussion_experiment_1';
save(filename)
