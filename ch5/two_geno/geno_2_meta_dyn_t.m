function[dx]=geno_2_meta_dyn_t(t,x);
% ode model for 2 genotype metapopulation dynamics

% ensure all variiables >= 0
x(1)=max(x(1),0);
x(2)=max(x(2),0);
x(3)=max(x(3),0);
x(4)=max(x(4),0);
% import parameters
global beta mu_M D_M psi1 u1 tu
% interpolate control - for cts control remove 'previous'
u=interp1(tu,u1,t,'previous');
% meta model when all means and variances >0
% f mean = 0 change in corresponding mean and variancce is zero
if x(1) >0 && x(2)>0
   dx=[(beta-mu_M-u)*x(1)-D_M*(x(1).^2+(1-psi1)*x(1)*x(2)+x(3)); (beta-mu_M)*x(2)-D_M*(x(2).^2+(psi1)*x(1)*x(2)+x(4));beta*x(1)-(mu_M+u)*(2*x(3)-x(1))-D_M*(4*x(3)^2/x(1)+4*x(1)*x(3)-3*x(3)-x(1)^2+(1-psi1)*(2*x(2)*x(3)-x(2)*x(1)));beta*x(2)-mu_M*(2*x(4)-x(2)) -D_M*(4*x(4)^2/x(2)+4*x(2)*x(4)-3*x(4)-x(2)^2+(psi1)*(2*x(1)*x(4)-x(2)*x(1)));];
elseif x(1)<=0 && x(2)<=0
   dx=[(beta-mu_M-u)*x(1)-D_M*(x(1).^2+(1-psi1)*x(1)*x(2)+x(3)); (beta-mu_M)*x(2)-D_M*(x(2).^2+(psi1)*x(1)*x(2)+x(4)); beta*x(1)-(mu_M+u)*(2*x(3)-x(1))-D_M*(4*x(1)*x(3)-3*x(3)-x(1)^2 +(1-psi1)*(2*x(2)*x(3)-x(2)*x(1))); beta*x(2)-mu_M*(2*x(4)-x(2))-D_M*(4*x(2)*x(4)-3*x(4)-x(2)^2+(psi1)*(2*x(1)*x(4)-x(2)*x(1)));];
elseif x(1)>0 && x(2)<=0
    dx=[(beta-mu_M-u)*x(1)-D_M*(x(1).^2+(1-psi1)*x(1)*x(2)+x(3)); (beta-mu_M)*x(2)-D_M*(x(2).^2+(psi1)*x(1)*x(2)+x(4)); beta*x(1)-(mu_M+u)*(2*x(3)-x(1))-D_M*(4*x(3)^2/x(1)+4*x(1)*x(3)-3*x(3)-x(1)^2 +(1-psi1)*(2*x(2)*x(3)-x(2)*x(1))); beta*x(2)-mu_M*(2*x(4)-x(2))-D_M*(4*x(2)*x(4)-3*x(4)-x(2)^2 +(psi1)*(2*x(1)*x(4)-x(2)*x(1)));];
else
     dx=[(beta-mu_M-u)*x(1)-D_M*(x(1).^2+(1-psi1)*x(1)*x(2)+x(3)); (beta-mu_M)*x(2)-D_M*(x(2).^2+(psi1)*x(1)*x(2)+x(4)); beta*x(1)-(mu_M+u)*(2*x(3)-x(1))-D_M*(4*x(1)*x(3)-3*x(3)-x(1)^2 +(1-psi1)*(2*x(2)*x(3)-x(2)*x(1))); beta*x(2)-mu_M*(2*x(4)-x(2))-D_M*(4*x(4)^2/x(2)+4*x(2)*x(4)-3*x(4)-x(2)^2 +(psi1)*(2*x(1)*x(4)-x(2)*x(1)));];
end
end
