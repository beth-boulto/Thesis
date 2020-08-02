function[dx]=geno_2_mean_dyn_t(t,x);
%ode system for two genotypes mean field model

% set means to be at least zero
x(1)=max(x(1),0);
x(2)=max(x(2),0);

% import parameters
global beta mu_M D_M psi1 u1 tu
% interpolate control - for cntinuous control remove previous
u=interp1(tu,u1,t,'previous');

   dx=[(beta-mu_M-u)*x(1)-D_M*(x(1).^2+(1-psi1)*x(1)*x(2));(beta-mu_M)*x(2)-D_M*(x(2).^2+(psi1)*x(1)*x(2))];

end
