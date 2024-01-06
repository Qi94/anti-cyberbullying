%{
Author: Y.Qin
Create: July.26.2022

to make sure state is in [0,1] after update
%}

function [output_State] = state_adjust(S,name)

adjust1 = (S>1); % S>1
adjust0 = (S<0); % S<0

if sum(adjust0 ~= 0)
    S(adjust0) = 0;
    disp([name,'<0']);
end

if sum(adjust1 ~= 0)
    S(adjust1) = 1;
    disp([name,'>1']);
end

output_State = S;

end

