% code to run basic analysis on three genotype model - simulation,eqm ,
%stability, response to perturbation
clear all
close all
% set parameters
global beta mu_M D_M psi1 psi2 psi3 u u2 tu gamma u1
beta = 1.5;
mu_M=1/(5*365);
D_M=1/20;
psi1=0.3;
psi2=0.1;
psi3=0.4;
% set ode model with defined function
dx=@(t,x)geno_3_dyn(t,x);
% set initial conditions [M_SS,M_SR,M_RR]
init=[0;5;0];
% set control intensity as zero
u=0;
u2=0;
% simulate for extened time -
[ts,xs]=ode45(dx,[0 150], init);
% approx eqm at end
eqm=xs(end,:);
% plot simulation to eqm
figure;
subplot(3,1,1)
plot(ts,xs(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('M_{ss}')
subplot(3,1,2)
plot(ts,xs(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('M_{sr}')
subplot(3,1,3)
plot(ts,xs(:,3),'g','LineWidth',4)
xlabel('time')
ylabel('M_{rr}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)



eqm=xs(end,:);

% calculate non-trivial jacobian with defined func
J_x= Jacob(eqm);
%  find eigenvalues
J_eig = eigs(J_x);
% find hermitian and eigenvalues for reactivity
H_x=1/2*(J_x+J_x');
H_eig=eigs(H_x);

% jacobian of single geno survivingeqms
J_ss=Jacob([(beta-mu_M)/D_M,0,0]);
J_ss_eig=eigs(J_ss);
H_ss=1/2*(J_ss+J_ss');
H_ss_eig=eigs(H_ss);


J_rr=Jacob([0,0,(beta-mu_M)/D_M]);

J_rr_eig=eigs(J_rr);
H_rr=1/2*(J_rr+J_rr');
H_rr_eig=eigs(H_rr);

% repeat for trivial eqm -from single genotype eqm 
J_0 = [beta-mu_M,0,0;0,0,0;0,0,0];
H_0=1/2*(J_0+J_0');
ei0=eigs(J_0);
Hei_0=eigs(H_0);

%% sustained treatment
% set initial conditiona t untreated eqm
init=eqm;
% apply constant treatment
u1=[1,1];
tu=[0, 50];
gamma=1;
u2=0;
mu_M=mu_M;
% run simulation with treatment for exteded time
dx=@(t,x)geno_3_dyn_t(t,x);
[ts2,xs2]=ode45(dx,[0 50],eqm);

% plot repsonse

figure;
subplot(3,1,1)
plot(ts2,xs2(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('M_{ss}')
subplot(3,1,2)
plot(ts2,xs2(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('M_{sr}')
subplot(3,1,3)
plot(ts2,xs2(:,3),'g','LineWidth',4)
xlabel('time')
ylabel('M_{rr}')
set(findall(gcf,'-property','FontSize'),'FontSize',20)

figure;
plot(ts2,(xs2(:,1)+xs2(:,2)+xs2(:,3)),'k','LineWidth',4)
xlabel('time')
ylabel('M_{total}')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
% calculate total mean
xt=xs2(:,1)+xs2(:,2)+xs2(:,3);
[tc,mc]=min(xt);

alp1=(eqm(1)+eqm(2)+eqm(3)-tc)/2;
pt=(2*xs2(:,3)+xs2(:,2))./(2*xt);

pm=min(pt);
pma=max(pt);

alpr=(pma-pm)/2;

%% post treat - follwng over treatment

% set init value to overtreated eqm
init=xs2(end,:);
% remove control
u1=[0,0];
tu=[0, 50];
gamma=0;
u2=0;
mu_M=mu_M;
% run simulation
dx=@(t,x)geno_3_dyn_t(t,x);
[ts2,xs2]=ode45(dx,[0 50],init);

% plot simulation resulst
figure;
subplot(3,1,1)
plot(ts2,xs2(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('M_{ss}')
subplot(3,1,2)
plot(ts2,xs2(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('M_{sr}')
subplot(3,1,3)
plot(ts2,xs2(:,3),'g','LineWidth',4)
xlabel('time')
ylabel('M_{rr}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

figure;
plot(ts2,(xs2(:,1)+xs2(:,2)+xs2(:,3)),'k','LineWidth',4)
xlabel('time')
ylabel('M_{total}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

xt=xs2(:,1)+xs2(:,2)+xs2(:,3);
[tc,mc]=min(xt);

alp1=(eqm(1)+eqm(2)+eqm(3)-tc)/2;
pt=(2*xs2(:,3)+xs2(:,2))./(2*xt);

pm=min(pt);
pma=max(pt);

alpr=(pma-pm)/2;
