%% code to numerically simulate mean field and variance inclusive models
%for a single genotype, plot results, determine stability and reactivity of
%equilibria, and plot phase plane of variance inclusive model.
clear all
close all
%set parameters
global beta  mu_L mu_M phi rho N D_M psi u tu
mu_M = 1.5-1/30;
D_M = 1/20;
beta=1.5;
mu_L=0.00014;
phi=0.005;
rho=0.504;
N=500;
psi=0.1;
% no treatment applied
u = [0,0];
tu = [0,1000];

% set initial conditions
Init_mean = 1;
Init_meta=[1;5];

% use mean field and metapop functions for ode models
dx=@(t,x)mean_field(t,x);
dx_v=@(t,x)metapop(t,x);
% run simulations
[t1,x1]=ode45(dx,[0 150],Init_mean);

[t2,x2]=ode45(dx_v,[0 10],Init_meta);


% comparison of mean plot
figure;
subplot(2,1,1)
plot(t1,x1,'r','Linewidth',4)
xlabel('time')
ylabel('mean sub-population size')
subplot(2,1,2)
plot(t2,x2(:,1),'r','Linewidth',4)
xlabel('time')
ylabel('mean sub-population size')
set(findall(gcf,'-property','FontSize'),'FontSize',20)

figure;
plot(t1,x1,'r','Linewidth',4)
xlabel('time')
ylabel('mean sub-population size')
set(findall(gcf,'-property','FontSize'),'FontSize',20)

%individual simulation plots
figure;
subplot(2,1,1)
plot(t2,x2(:,1),'r','Linewidth',4)
xlabel('time')
ylabel('mean sub-population size')
subplot(2,1,2)
plot(t2,x2(:,2),'b','Linewidth',4)
xlabel('time')
ylabel('variance in sub-population size')
set(findall(gcf,'-property','FontSize'),'FontSize',20)

%% phase plane for variance model
M=0:30;
% find variance as function of mean from mean eq
V1 = ((beta-mu_M)*M-D_M*M.^2)/D_M;

% find variance as function of mean from var eq

A = -4*D_M./M;
B = -2*mu_M+3*D_M-4*M*D_M;
C =beta*M+mu_M*M+D_M*M.^2;

V2 = (-B+(B.^2-4.*A.*C).^(1/2))./(2*A);
V3 =(-B-(B.^2-4.*A.*C).^(1/2))./(2*A);
% set multiple initial cnditions for trajectories
M_in=[ones(1,5),10*ones(1,5),15*ones(1,5),25*ones(1,5)];
V_in = [10,25,50,75,100];
V_in=horzcat(V_in,V_in,V_in,V_in);
% set up mean and variance as meshgrid
[x1,y1]=meshgrid(0.00001:01:30,0.01:1:30);
u1=zeros(size(x1));
v1=zeros(size(x1));
% calculate derivatives for each combo of M and V in meshgrid
for ii=1:numel(x1)
     u1(ii) = beta*x1(ii)-(mu_M)*x1(ii)-D_M*(x1(ii).^2+y1(ii));
  v1(ii)=beta*x1(ii)-mu_M*(2*y1(ii)-x1(ii))-D_M*(4*y1(ii)^2/x1(ii)-3*y1(ii)-x1(ii)^2+4*x1(ii)*y1(ii));
  u1(ii)=u1(ii)./((u1(ii)^2+v1(ii)^2)^(1/2));
   v1(ii)=v1(ii)./((u1(ii)^2+v1(ii)^2)^(1/2));
end
% plot nullclines
figure;

plot(M,V1,'r',M,V2,'b',M,V3,'b',M,M,'k','LineWidth',4)
xlabel('mean')
ylabel('variance')
hold on
%plot derivatives at all ppoints
quiver(x1,y1,u1,v1)
% plot trajectories from diff initial conditions
for ii=1:20
    [t3,x3]=ode45(dx_v,[0 50], [M_in(ii);V_in(ii)]);
    plot(x3(:,1),x3(:,2),'g','LineWidth',4)
end
set(findall(gcf,'-property','FontSize'),'FontSize',20)


%% Jacobians

% mean field jacobian - non trivial state

J = [beta-mu_M-2*D_M*(beta-mu_M)/D_M];
%eigenvalues
ei=eigs(J);
%resilience
Res_M = -real(ei(1));
%reactivity
H =1/2*J+1/2*transpose(J);
Hei=eigs(H);
Hei = real(Hei);
Hei = sort(Hei,'descend');
Reac_M =Hei(1);
% trivial state
J_0 = [beta-mu_M];
ei_0=eigs(J_0);

% variance inclusive

J_v = [beta-mu_M-2*D_M*x2(end,1),-D_M;
    beta+mu_M-D_M*(-4*x2(end,2).^2/(x2(end,1).^2)-2*x2(end,1)+4*x2(end,2)), -2*mu_M-D_M*(8*x2(end,2)/x2(end,1)-3+4*x2(end,1))];
ei_v=eigs(J_v);

J_0 =  [beta-mu_M-2*D_M*0,-D_M;
    beta+mu_M-D_M*(-2*0+4*0), -2*mu_M-D_M*(-3+4*0)];
ei_0v=eigs(J_0);
H_0=1/2*(J_0+J_0');
Hei_0=eigs(H_0);

ei_v=real(ei_v);
ei_v=sort(ei_v,'descend');
%resilience
Res_v = -(ei_v(1));

H_v=1/2*J_v+1/2*transpose(J_v);
Hei_v=eigs(H_v);
Reac_v =Hei_v(1);

%% Rsponse to perturbation

% perturb from equilibrium
Init_m = x1(end) - 5;

Init_v = [x2(end,1)-5;x2(end,2)];

% simulate return to eqm
[t1,x1]=ode45(dx,[0 2], Init_m);
[t2,x2]=ode45(dx_v,[0 2], Init_v);

% plot behaviour
figure;
plot(t1,x1,'r','Linewidth',4)
xlabel('time')
ylabel('mean sub-population size')
set(findall(gcf,'-property','FontSize'),'FontSize',20)

figure;
subplot(2,1,1)
plot(t2,x2(:,1),'r','Linewidth',4)
xlabel('time')
ylabel('mean sub-population size')
subplot(2,1,2)
plot(t2,x2(:,2),'b','Linewidth',4)
xlabel('time')
ylabel('variance in sub-population size')
set(findall(gcf,'-property','FontSize'),'FontSize',20)

%% eigenvalue sof rivial eqm variance model

J_v = [beta-mu_M-2*D_M*0, -D_M;
    beta+mu_M-D_M*(-2*0), -mu_M-3*D_M];
ei_v=eigs(J_v);

ei_v=real(ei_v);
ei_v=sort(ei_v,'descend');
%resilience
Res_v = -(ei_v(1));

H_v=1/2*J_v+1/2*transpose(J_v);
Hei_v=eigs(H_v);
Hei_v = real(Hei_v);
Hei_v = sort(Hei_v,'descend');
Reac_v =Hei_v(1);

J_v2 = [beta-mu_M-2*D_M*x2(end,1), -D_M;
    beta+mu_M-D_M*(-4*x2(end,2)^2/(x2(end,1).^2)+4*x2(end,2)-2*x2(end,1)),-2*mu_M-D_M*(8*x2(end,2)/x2(end,1)-3+4*x2(end,1))];

ei_v2=eigs(J_v2);
H_v2 = 1/2*(J_v2+transpose(J_v2));
Hei_v2=eigs(H_v2);


J_v = [beta-mu_M-2*D_M*x2(end,1),-D_M;
    beta+mu_M-D_M*(-4*x2(end,2).^2/(x2(end,1).^2)-2*x2(end,1)+4*x2(end,2)), -2*mu_M-D_M*(8*x2(end,2)/x2(end,1)-3+4*x2(end,1))];
H_v=1/2*J_v+1/2*transpose(J_v);
Hei_v=eigs(H_v);