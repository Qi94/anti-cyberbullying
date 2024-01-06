function [gu] = SCR_g_b_normal(u)
%a function that is neither conver nor concave 
% f(u)=2/(1+e^(1-U))

%gu = 2 ./ (1.2 + exp(1-u));
gu = 0.14.*sqrt(u);
end