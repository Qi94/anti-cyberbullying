function [vu] = SCR_v_normal(u)
%a function that is neither conver nor concave 
% f(u)=2/(1+e^(1-U))

%vu = 2 ./ (1 + exp(1-u));
vu = 0.12 ./ (1 + exp(1-u));
end