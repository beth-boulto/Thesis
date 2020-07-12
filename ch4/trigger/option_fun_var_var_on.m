function [y,isterm,dir ] = option_fun_var_var_on(t,x)
% option function to start treatment on variance inclusive model when mean
% variance get high enough. event marked by change in y value from -1 to
% 1
%set global variables
global phi alpha mu_M p r u1 tu MT VT T_start
dx_m=@(t,x)metapop(t,x);


% if treatment has not just ended
if t>T_start+0.0001
    % if v high enough set event occurance to yes ( y=1)
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
