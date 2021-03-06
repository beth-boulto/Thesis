function [y,isterm,dir ] = option_fun_ish_both_off(t,x)
% event function for isham model occurs if both mean and variance below
% set value - event is shown by change in y value from y=-1 to y=1.
%set global variables
global beta mu_M D_M psi1 MT  T_start
dx_m=@(t,x)isham(t,x);
change_crit=1;
m=dx_m(t,x);
m1=abs(m(1));
m2=abs(m(2));
% if treatment has not just startedd
if t>T_start+0.0001
    %check if mean and variance below value
if m1<change_crit && m2<change_crit
    % if yes set y=1
y = 1;
else 
    %otherwise y=-1
y=-1  ;
end
else
%y value set to be negative    
    y = -1;
end
isterm=1;
dir=0;

end
