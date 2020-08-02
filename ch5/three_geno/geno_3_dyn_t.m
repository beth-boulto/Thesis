function[dx]=geno_3_dyn_t(t,x);

% three genotype ode model function

% ensure all means >=0
x(1)=max(x(1),0);
x(2)=max(x(2),0);
x(3)=max(x(3),0);
% import parameters
global beta mu_M D_M psi1 psi2 psi3 tu u1 gamma
% interpolate for control - for sequential set to 'previous' for cts remove
% option
u = interp1(tu,u1,t,'previous');
% set control on mixed geno
u2 = gamma*u;
% if all eradicated no change
if x(1)+x(2)+x(3)==0
    dx=[0;0;0];
else
    A=(2*x(1).*x(2)+x(1).*x(2)+x(2).*x(3)+0.5*x(2).^2);
     dx=zeros(3,1);
   dx(1)= beta*(x(1).^2+x(1).*x(2)+0.25*x(2).^2)./(x(1)+x(2)+x(3))-(mu_M+u)*x(1)-D_M*(x(1).^2+psi1*x(1).*x(2)+psi2*x(1).*x(3));
   dx(2) = beta*A./(x(1)+x(2)+x(3))-(mu_M+u2)*x(2)-D_M*(x(2).^2+(1-psi1)*x(1).*x(2)+psi3*x(2).*x(3));
   dx(3) = beta*(x(3).^2+x(3).*x(2)+0.25*x(2).^2)./(x(1)+x(2)+x(3)) -mu_M*x(3)-D_M*(x(3).^2+(1-psi2)*x(1).*x(3)+(1-psi3)*x(3).*x(2));
   
end
end
