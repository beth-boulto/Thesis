function[J]=Jacob(x)
% sets up jacobian of three genotype model for non-trivial eqm
global beta mu_M D_M psi1 psi2 psi3
A=x(1)+x(2)+x(3);
J(1,1)=((A*(2*x(1)+x(2))-(x(1)^2+x(1)*x(2)+0.25*x(2)^2))/(A^2))-mu_M-D_M*(2*x(1)+psi1*x(2)+psi2*x(3));
J(1,2)=(((x(1)+0.5*x(2))*A-x(1)^2+x(1)*x(2)+0.25*x(2)^2)/(A^2))-D_M*(psi1*x(1));
J(1,3)=(-(x(1)^2+x(1)*x(2)+0.25*x(2).^2)/(A^2))-D_M*(psi2*x(1));
J(2,1)=((A*(x(3)+x(2))-(x(1)*x(2)+x(1)*x(3)+0.5*x(2)^2+x(2)*x(3)))/(A^2))-D_M*(1-psi1)*x(2);
J(2,2)=((A*(A)-(x(1)*x(2)+x(1)*x(2)+0.5*x(2)^2+x(3)*x(2)))/(A^2))-mu_M-D_M*(2*x(2)+(1-psi1)*x(1)+psi3*x(3));
J(3,1)=((A*(x(1)+x(2))-(x(1)*x(2)+x(1)*x(3)+0.5*x(2)^2+x(2)*x(3)))/(A^2))-D_M*(psi3)*x(2);
J(3,1)=(-(x(3)^2+x(3)*x(2)+0.25*x(2).^2)/(A^2))-D_M*(1-psi2)*x(3);
J(3,2)=(((x(3)+0.5*x(2))*A-x(3)^2+x(3)*x(2)+0.25*x(2)^2)/(A^2)) -D_M*((1-psi3)*x(3));
J(3,3)=((A*(2*x(3)+x(2))-(x(1)^2+x(3)*x(2)+0.25*x(2)^2))/(A^2))-mu_M-D_M*(2*x(3)+(1-psi3)*x(2)+(1-psi2)*x(1));
J=beta*J;
end
