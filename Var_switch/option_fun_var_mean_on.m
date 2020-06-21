function [y,isterm,dir ] = option_fun_var_mean_on(t,x)
% option function to start treatment on variance inclusive model when mean
% gets high enough. event marked by change in y value from -1 to
% 1
%set global variables
global beta mu_M D_M psi1 MT  T_start
dx_m=@(t,x)metapop(t,x);


% if treatment has not just ended
if t>T_start+0.0001
    % check if mean is high enough. if yes y=1, if no y==-1
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
