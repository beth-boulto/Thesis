function[J]=jacob_non_triv(x)
% code to calculte the jacobian of two geno variance model at non-trivial eqm
% import parameters
global beta mu_M D_M psi1 
J=zeros(4,4);
% when M_U>0
if x(1)>0
    J(3,1) = beta+mu_M-D_M*(-4*x(3)^2/(x(1)^2)+4*x(3)-2*x(1)+(1-psi1)*(-x(2)));
    J(3,3)=-2*mu_M-D_M*(8*x(3)/x(1)+4*x(1)-3+(1-psi1)*(x(2)));
else %M_U<=0
    J(3,1)=beta+mu_M-D_M*(4*x(3)-2*x(1)+(1-psi1)*(-x(2)));
    J(3,3)=-2*mu_M-D_M*(4*x(1)-3+(1-psi1)*(x(2)));
end
if x(2)>0 %M_T>0
    J(4,2)=beta+mu_M-D_M*(-4*x(4)^2/(x(2)^2)+4*x(4)-2*x(2)+(psi1)*(-x(1)));
    J(4,4)=-2*mu_M-D_M*(8*x(4)/x(2)+4*x(2)-3+(psi1)*(x(1)));

else %M_T<=0
     J(4,2)=beta+mu_M-D_M*(4*x(4)-2*x(2)+(psi1)*(-x(1)));
    J(4,4)=-2*mu_M-D_M*(4*x(2)-3+(psi1)*(x(1)));
end
% where M_U & M_T dont matter
J(1,1) = beta-mu_M-D_M*(2*x(1)+(1-psi1)*x(2));
J(1,2)=-D_M*(1-psi1)*x(1);
J(1,3) = -D_M;
J(1,4)=0;
J(2,1)=-D_M*psi1*x(2);
J(2,2)=beta-mu_M-D_M*(2*x(2)+(psi1)*x(1));
J(2,3)=0;
J(2,4)=-D_M;

J(3,2)=-D_M*((1-psi1)*(2*x(3)-x(1)));

J(3,4)=0;
J(4,1)=-D_M*((psi1)*(2*x(4)-x(2)));

J(4,3)=0;



end