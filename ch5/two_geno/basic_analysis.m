%% code for basic analysis of two genotype models -simulation, equilibrium, 
%jacobian stability and reactivity
clear all
close all
% set parameters
global beta mu_M D_M psi1  u
beta=1.5;
mu_M= 1/(5*365);
D_M=1/20;
u=0;
psi1=0.9;
%mu_M=1.6;
 % set ode models for mean field model and metapop model using predefined
 % functions
dx_m=@(t,x)geno_2_mean_dyn(t,x);

dx_v=@(t,x)geno_2_meta_dyn(t,x);
% set initial conditions
init_m=[25;5];
init_v=[1;1;1;1];
% run for extended time to determine stability
[tm,xm]=ode45(dx_m, [0 15],init_m);
[tv,xv]=ode45(dx_v, [0 500],init_v);
% plot simulation to eqm
figure;
subplot(2,1,1)
plot(tm,xm(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('M_{s}')
subplot(2,1,2)
plot(tm,xm(:,2),'r','LineWidth',4)
xlabel('time')
ylabel('M_{r}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)


figure;
subplot(2,2,1)
plot(tv,xv(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('M_{s}')
subplot(2,2,3)
plot(tv,xv(:,2),'r','LineWidth',4)
xlabel('time')
ylabel('M_{r}')
subplot(2,2,2)
plot(tv,xv(:,3),'--r','LineWidth',4)
xlabel('time')
ylabel('V_s')
subplot(2,2,4)
plot(tv,xv(:,4),'--b','LineWidth',4)
xlabel('time')
ylabel('V_r')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

%% JAcobian stability analysis
% set jacobian matrix via predefined function and aproximate eqm
J_v=jacob_non_triv(xv(end,:));
% determine eigenvalues
ei_v=eigs(J_v);
% calculate hermitian matrix
H_v=1/2*(J_v+J_v');
% find eigenvalues
Hei_v=eigs(H_v);

% set jacobiaan when M=V=0 eqm
J_0 = [beta-mu_M,0,-D_M,0;0,beta-mu_M,0,-D_M;beta+mu_M,0,-D_M*(-3),0;
    0,beta+mu_M,0,-D_M*(-3)];
% find eigenvalues
ei_0=eigs(J_0);

% find hermitian and eigenvalues
H_0=1/2*(J_0+J_0');
Hei_0=eigs(H_0);


%% sustained treatment - repsosne
% set treatment as consatn

u=1;

% set initial conditions as approx eqms
init_m=xm(end,:);
init_v=xv(end,:);

% simulate for extended period
[tm,xm]=ode45(dx_m,[0 30],init_m);

[tv,xv]=ode45(dx_v,[0 30],init_v);

% calculate total mean for mf model
xmt=xm(:,1)+xm(:,2);
[mm,mt]=min(xmt);

% calculate total mean and variance for v model
xvm=xv(:,1)+xv(:,2);
xvv=xv(:,3)+xv(:,4);

% find when minima are reached
[mvm,mvt]=min(xvm);
[vm,vt]=min(xvv);
 
% plot responses to control
figure;
subplot(2,1,1)
plot(tm,xm(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('M_{s}')
subplot(2,1,2)
plot(tm,xm(:,2),'r','LineWidth',4)
xlabel('time')
ylabel('M_{r}')
set(findall(gcf,'-property','FontSize'),'FontSize',20)

figure;
plot(tm,xm(:,1)+xm(:,2),'k','LineWidth',4)
xlabel('time')
ylabel('M_{total}')
set(findall(gcf,'-property','FontSize'),'FontSize',20)

figure;
subplot(2,2,1)
plot(tv,xv(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('M_{s}')
subplot(2,2,2)
plot(tv,xv(:,2),'--r','LineWidth',4)
xlabel('time')
ylabel('M_{r}')
subplot(2,2,3)
plot(tv,xv(:,3),'b','LineWidth',4)
xlabel('time')
ylabel('V_s')
subplot(2,2,4)
plot(tv,xv(:,4),'--b','LineWidth',4)
xlabel('time')
ylabel('V_r')
set(findall(gcf,'-property','FontSize'),'FontSize',20)

figure;
subplot(2,1,1)
plot(tv,xv(:,1)+xv(:,2),'k','LineWidth',4)
xlabel('time')
ylabel('M_{total}')

subplot(2,1,2)
plot(tv,xv(:,3)+xv(:,4),'--k','LineWidth',4)
xlabel('time')
ylabel('V_{total}')
set(findall(gcf,'-property','FontSize'),'FontSize',20)

%% post treatment - response to overtreatment
% set control to constant 0
u=0;

% set initial condiitons to end conditionsof sustained treatment
init_m=xm(end,:);
init_v=xv(end,:);

% simulate for extended period
[tm,xm]=ode45(dx_m,[0 30],init_m);

[tv,xv]=ode45(dx_v,[0 30],init_v);

% calculate total variables
xmt=xm(:,1)+xm(:,2);
[mm,mt]=min(xmt);

xvm=xv(:,1)+xv(:,2);
xvv=xv(:,3)+xv(:,4);

[mvm,mvt]=min(xvm);
[vm,vt]=min(xvv);
 
% plot simulations
figure;
subplot(2,1,1)
plot(tm,xm(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('M_{s}')
subplot(2,1,2)
plot(tm,xm(:,2),'r','LineWidth',4)
xlabel('time')
ylabel('M_{r}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

figure;
plot(tm,xm(:,1)+xm(:,2),'k','LineWidth',4)
xlabel('time')
ylabel('M_{total}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

figure;
subplot(2,2,1)
plot(tv,xv(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('M_{s}')
subplot(2,2,2)
plot(tv,xv(:,2),'--r','LineWidth',4)
xlabel('time')
ylabel('M_{r}')
subplot(2,2,3)
plot(tv,xv(:,3),'b','LineWidth',4)
xlabel('time')
ylabel('V_s')
subplot(2,2,4)
plot(tv,xv(:,4),'--b','LineWidth',4)
xlabel('time')
ylabel('V_r')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

figure;
subplot(2,1,1)
plot(tv,xv(:,1)+xv(:,2),'k','LineWidth',4)
xlabel('time')
ylabel('M_{total}')

subplot(2,1,2)
plot(tv,xv(:,3)+xv(:,4),'--k','LineWidth',4)
xlabel('time')
ylabel('V_{total}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

