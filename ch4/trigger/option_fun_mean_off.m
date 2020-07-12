function [y,isterm,dir ] = option_fun_mean_off(t,x)
% event function to stop treatment on mean field model when mean is low
% enough. event occurs when y changes from 1 to -1.
%set global variables
global beta mu_M D_M psi1 MT  T_start
dx_m=@(t,x)mean_field(t,x);
change_crit=1;
m=dx_m(t,x);
m=abs(m);
% if treatment has not juststarted
if t>T_start+0.0001
    % check if mean low enough
if m<=change_crit
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
