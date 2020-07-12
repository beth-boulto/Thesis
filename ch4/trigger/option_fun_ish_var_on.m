function [y,isterm,dir ] = option_fun_ish_var_on(t,x)
% event function for isham model occurs if variance above
% set value - event is shown by change in y value from y=-1 to y=1.
%set global variables
global phi alpha mu_M p r u1 tu MT VT  T_start
dx_m=@(t,x)isham(t,x);

% if treatment has not just ended
if t>T_start+0.0001
    % if variance too high simulation ends
if x(2)>VT 
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
