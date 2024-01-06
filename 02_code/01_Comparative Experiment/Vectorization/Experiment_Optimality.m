load plot-twitter10000.mat
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

filename = 'optimalResult-twitter10000.mat';
save(filename)