function [y,isterm,dir ] = option_fun_var_mean_off(t,x)
% option function to stop treatment on variance inclusive model when mean
% get low enough. event marked by change in y value from -1 to
% 1
%set global variables
global beta mu_M D_M psi1 MT  T_start
dx_m=@(t,x)metapop(t,x);
change_crit=1;
m=dx_m(t,x);
m=abs(m(1));
% if treatment has not just ended
if t>T_start+0.0001
    % check if mean is low enough - if yes set y=1. else y=-1
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
