function [y,isterm,dir ] = option_fun_ish_prop_on(t,x)
% event function for isham model occurs if both proportion of hosts
% exceeding set burden exceeds
% set value - event is shown by change in y value from y=-1 to y=1.
%set global variables
global phi alpha mu_M p r prop  T_start C_num
dx_m=@(t,x)isham(t,x);
% set up approximate distribution
p1=x(1)/x(2);
r1=x(1)^2/(x(2)-x(1));

pd1=makedist('NegativeBinomial','p',p1,'R',r1);
p2=cdf(pd1,C_num,'upper');

% if treatment has not just ended
if t>T_start+0.0001
    % check if p2 is too high
if p2>prop
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
