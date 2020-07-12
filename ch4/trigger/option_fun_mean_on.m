function [y,isterm,dir ] = option_fun_mean_on(t,x)
% event function to start treatment on mean field model when mean is high
% enough. event occurs when y changes from 1 to -1.
%set global variables
global beta mu_M D_M psi1 MT  T_start
dx_m=@(t,x)mean_field(t,x);


% if treatment has not just ended
if t>T_start+0.0001
    % check if mean high enough
if x(1)>MT 
y = 1;
else 
y=-1  ;
end
else
%y value set to be negative    
    y = -1;
end
isterm=1;
dir=0;

end
