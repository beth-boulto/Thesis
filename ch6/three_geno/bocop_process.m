% bocop resulst plotting code for three geno model - imports data frombocop
% optimsiation and plots results and calculates total cost
clear all
close all
%set parameters
global beta mu_M D_M psi1 psi2 psi3 u1 tu gamma
beta = 1.5;
mu_M= 1/(5*365);
D_M=1/20;
psi1=0.3;
psi2=0.1;
psi3=0.4;
gamma=0;
% set ode model usig defined function
dx=@(t,x)geno_3_dyn_t(t,x);
% set initial conditions same as bocop i.c's
init=[23.2872,   28.9256,    4.4197];
%import bocop data
bocop1 = readtable('three_geno_mod.xlsx');
us1 = bocop1{1:end,{'paramset1'}};
us2 = bocop1{1:end,{'paramset2'}};
us3 = bocop1{1:end,{'paramset3'}};
us4 = bocop1{1:end,{'paramset4'}};
% choose control
u1=us4;
% set time vector for control
m=length(u1);
tu=0:30/(m-1):30;

% simulate system with chosen control
[tm,xm]=ode45(dx,[0 30],init);
% plotcontrol and resposnse 
figure;
subplot(2,1,1)
plot(tu,u1,'k','LineWidth',4)
xlabel('Time')
ylabel('control intensity')



subplot(2,1,2)
plot(tm,xm(:,1)+xm(:,2)+xm(:,3),'k','LineWidth',4)
xlabel('Time')
ylabel('total mean population')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

figure;
subplot(3,1,1)
plot(tm,xm(:,1),'r','LineWidth',4)
xlabel('Time')
ylabel('M_{ss}')
subplot(3,1,2)
plot(tm,xm(:,2),'b','LineWidth',4)
xlabel('Time')
ylabel('M_{sr}')
subplot(3,1,3)
plot(tm,xm(:,3),'g','LineWidth',4)
xlabel('Time')
ylabel('M_{rr}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
