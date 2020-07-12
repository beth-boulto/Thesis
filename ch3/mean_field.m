function[dx]=mean_field(t,x)
%% MEAN FIELD ODE MODEL
% import parameters
global beta mu_M D_M tu u
% interpolate for treatment intensity
u1=interp1(tu,u,t);
% set dM/dt
dx=beta*x(1)-(mu_M+u1)*x(1)-D_M*x(1)^2;

end