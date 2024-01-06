function [fu] = SCR_f_b_normal(u)
%a function that is neither conver nor concave 
% f(u)=3/(1+e^(1-U))

%fu = 3 ./ (1 + exp(1-u));
fu = 0.15.*sqrt(u);
end
