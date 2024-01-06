function [gu] = SCR_g_c_normal(u)
%a function that is neither conver nor concave 
% f(u)=2/(1+e^(1-U))
%gu = 2 ./ (1.1 + exp(1-u));
gu = 0.18.*sqrt(u);
end