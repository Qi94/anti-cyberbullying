function [fu] = SCR_f_a_normal(u)
%a function that is neither conver nor concave 
% f(u)=2/(1+e^(1-U))
%fu = 2 ./ (1 + exp(1-u));
fu = 0.13.*sqrt(u);
end
