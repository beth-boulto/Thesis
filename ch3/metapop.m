function[dx]=metapop(t,x)
%% reusable function for metapopulation ode model

% import global parameters
global beta mu_M D_M tu u
% ensure mean and variance do not decrease below zero
x(1)=max(x(1),0);
x(2)=max(x(2),0);
% interpolate to find current treatment intensity
u1=interp1(tu,u,t);
% when variance is greater than zero model
if x(1)>0 && x(2)>=0
dx=[beta*x(1)-(mu_M+u1)*x(1)-D_M*(x(2)+x(1)^2);beta*x(1)-(mu_M+u1)*(2*x(2)-x(1))-D_M*(4*x(2)^2/x(1)-3*x(2)-x(1)^2+4*x(1)*x(2))];
elseif x(1)>0 && x(2)<0
  % variance less than zero model - MEAN FILED MODEL
    dx=[beta*x(1)-(mu_M+u1)*x(1)-D_M*(x(1)^2);beta*x(1)-(mu_M+u1)*(-x(1))-D_M*(x(1)^2)];
else
    %when mean and variance both less than or equal to zero no chNGE
    dx=[0;0];
end

end
