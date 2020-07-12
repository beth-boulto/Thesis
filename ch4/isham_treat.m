function[dx]=isham_treat(t,x)
%% isham ode model function for use in treated simulations
% import parameters
global phi alpha mu_M p r u1 tu
% interpolate for current control
u=interp1(tu,u1,t);
% calculate neg bin pgf plus derivatives
h1=h_neg_bin(1);
% mean and variance derivatives
dx=zeros(2,1);
dx(1)=phi * (h1(2))  - alpha* x(2) - (mu_M+u)*x(1);
dx(2)=phi*(h1(3)+ h1(2))+ (mu_M+u)*x(1)-2*(mu_M+u)*x(2);
 end