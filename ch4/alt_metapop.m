function[dx]=alt_metapop(t,x)
% alternate metapopulation ode for controll application

% import parameters
global beta mu_M D_M tu u
x(1)=max(x(1),0);
x(2)=max(x(2),0);
% find current control intensity by previous value
u1=interp1(tu,u,t,'previous');
% variance model with V>=0
if x(1)>0 && x(2)>=0
dx=[beta*x(1)-(mu_M+u1)*x(1)-D_M*(x(2)+x(1)^2);beta*x(1)-(mu_M+u1)*(2*x(2)-x(1))-D_M*(4*x(2)^2/x(1)-3*x(2)-x(1)^2+4*x(1)*x(2))];
elseif x(1)>0 && x(2)<0
  % mean field model when V<=0
    dx=[beta*x(1)-(mu_M+u1)*x(1)-D_M*(x(1)^2);beta*x(1)-(mu_M+u1)*(-x(1))-D_M*(x(1)^2)];
else
    % no change when both M & V>=0
    dx=[0;0];
end

end
