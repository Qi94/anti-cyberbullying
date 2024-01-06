clear;
clc;
%load 'optimalResult-Youtube100.mat'
load 'optimalResult-Youtube10000.mat' 

%% ===============Comparative Experiment===============

% ---initiate array to store the comparative experiments---
experiment_times = 100;
J_Y = zeros(1, experiment_times);

% ---keep the optimality result---
R_x_optimality = r_x;

% ---calculate anti-cyberbullying loss---
iter = 1;
while true
    % --preliminary preparation--
    % report the number of experiment
    disp(['progressing ',num2str(iter),'th experiment...']);

    % ---randomize the R_Y when cyberbullying insiste on the R_X^*---
    t_1_y = randi([0 n],1,1);
    t_2_y = randi([0 n],1,1);

    % ensure the time series for the random time point
    if t_1_y > t_2_y
        y = t_2_y;
        t_2_y = t_1_y;
        t_1_y = y;
    end

    % initiate benefit and loss parameters for strategy profile
    J_y_iter = 0;

    % ---Revert back to original OSN state---
    E_C = zeros(N, n+1);
    E_A = zeros(N, n+1);
    E_C(:,1) = 0.3*ones(N,1);
    E_A(:,1) = 0.3*ones(N,1);

    % ---initiate new array for random strategy---
    r_y_random = zeros(1, n+1);
    
    % generate radom strategy
    for t = 1:n
        if t < t_1_y
            r_y_random(1,t) = r_y_max;
            
        elseif (t >= t_1_y) && (t <= t_2_y)
            r_y_random(1, t) = (r_y_max)*((t_2_y - t)/(t_2_y - t_1_y));
            
        elseif t > t_2_y
            r_y_random(1,t) = 0;

        end
    end

    % calculate loss for anti-cyberbullying
    %1.random r_y
    %2.forward propagate -> E_C
    %3.r_y & E_C -> J_y
    for t = 1:n
        % --calculate the state evolvement during each iteration--

        % calculate each node separately for simplicity
        dc = (SCR_f_b_normal(r_x(t)) + beta(2) * a * E_C(:,t) + SCR_m_normal(r_y_random(t)) + zeta * a * E_A(:,t)) .* ...
            (1 - E_C(:,t) - E_A(:,t)) + (SCR_f_a_normal(r_x(t)) + beta(1) * a * E_C(:,t)) .* E_A(:,t) - ...
            (SCR_g_c_normal(r_y_random(t)) + gamma(2) * a * E_A(:,t)) .* E_C(:,t);
        
        da = (SCR_g_b_normal(r_y_random(t)) + gamma(1) * a * E_A(:,t) + SCR_v_normal(r_x(t)) + eta * a * E_C(:,t)) .* ...
            (1 - E_C(:,t) - E_A(:,t)) + (SCR_g_c_normal(r_y_random(t)) + gamma(2) * a * E_A(:,t)) .* E_C(:,t) - ...
            (SCR_f_a_normal(r_x(t)) + beta(1) * a * E_C(:,t)) .* E_A(:,t);
        
        % update
        E_C(:,t+1) = E_C(:,t) + dt * dc;
        E_A(:,t+1) = E_A(:,t) + dt * da;
        
        % make sure E_A and E_C is in [0,1]
        E_C(:,t+1) = state_adjust(E_C(:,t+1), 'c');
        E_A(:,t+1) = state_adjust(E_A(:,t+1), 'a');  
    
        term_J = 0;
        for i = 1:N
            term_J = term_J + w(2)*E_C(i,t);
        end
        J_y_iter = J_y_iter + term_J*dt + r_y_random(1,t)*dt;
    end

    % --save the total benefit for the random anti-cyberbulluing stategy
    J_Y(1, iter) = J_y_iter;
    
    % -- start a new experimenthaode 
    iter = iter+1;

    % --check if the experiment times are enough
    if iter > experiment_times
    disp(['--algorithm terminated at ',num2str(iter),'th experiment--']);
    break;
    end
end


filename = 'Experiment-anti-cyberbullying-Youtube100.mat';
save(filename)
