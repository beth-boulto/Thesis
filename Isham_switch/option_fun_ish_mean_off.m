function [y,isterm,dir ] = option_fun_ish_mean_off(t,x)
% event function for isham model occurs if mean below
% set value - event is shown by change in y value from y=-1 to y=1.
%set global variables
global beta mu_M D_M psi1 MT  T_start
dx_m=@(t,x)isham(t,x);
change_crit=1;
m=dx_m(t,x);
m1=abs(m(1));
m2=abs(m(2));
% if treatment has not just started
if t>T_start+0.0001
if m1<change_crit 
    % check if mean less than set value
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
