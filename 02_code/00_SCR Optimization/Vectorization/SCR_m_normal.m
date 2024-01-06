function [mu] = SCR_m_normal(u)
%a function that is neither conver nor concave 
% f(u)=2/(1+e^(1-U))

%mu = 2 ./ (1 + exp(1-u));
mu = 0.13 ./ (1 + exp(1-u));
end
