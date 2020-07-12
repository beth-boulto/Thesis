%% code to import the data from the bocop optimisation of the Isham model 
% and plot the results 
clear all
close all
global phi alpha mu_M p r u1 tu
% set parameter values
phi = 52/12;
mu_H = 0;
mu_M= 1/12;
u=20;
alpha=0.002;
% set mean clump size
M = 0.5;
% set dispersion parameter

    r=0.5;    
% fin variance in clump sizes

   
V = (M.^2+r*M)/r;
   %find p parameter 
  prob=M/V;
  p=1-prob;
  
  %import bocop data
bocop1 = readtable('isham_halving.xlsx');
u_m = bocop1{1:end,{'mean'}};
u_vm=bocop1{1:end,{'meanvar'}};
u_v=bocop1{1:end,{'var'}};
m=length(u_vm);
tu=0:10/(m-1):10;
dx= @(t,x)isham_treat(t,x);

%% mean control
% set control to bocop import
u1=u_m;
% simulate with isham ode model
[tm,xm]=ode45(dx,[0 10],[25,45]);
% calculate cost via cost function
Y = 0.00754*(xm(:,1)-1.98).^2;
costm= trapz(tu,u1);
costm=costm+trapz(tm,Y);
costm2=trapz(tu,u1);
%plot control and response
figure;
subplot(3,1,1)
plot(tu,u1,'k','LineWidth',4)
xlabel('time')
ylabel('control intensity')
subplot(3,1,2)
plot(tm,xm(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('mean population size')
subplot(3,1,3)
plot(tm,xm(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('population size variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)


%% variance control
% set control to bocop imprt
u1=u_v;
% simulate on isham odemodel
[tv,xv]=ode45(dx,[0 10],[25,45]);
%calculate cost function
Y =0.00238*(xv(:,2)-3.49).^2;
costv= trapz(tu,u1);
costv=costv+trapz(tv,Y);
costv2=trapz(tu,u1);
%plot control and resposne
figure;
subplot(3,1,1)
plot(tu,u1,'k','LineWidth',4)
xlabel('time')
ylabel('control intensity')
subplot(3,1,2)
plot(tv,xv(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('mean population size')
subplot(3,1,3)
plot(tv,xv(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('population size variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)


%% both control

u1=u_vm;

[tb,xb]=ode45(dx,[0 10],[25,45]);
Y = 0.00377*(xb(:,1)-1.98).^2+0.00119*(xb(:,2)-3.49).^2;
costb= trapz(tu,u1);
costb=costb+trapz(tb,Y);
costb2=trapz(tu,u1);
figure;
subplot(3,1,1)
plot(tu,u1,'k','LineWidth',4)
xlabel('time')
ylabel('control intensity')
subplot(3,1,2)
plot(tb,xb(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('mean population size')
subplot(3,1,3)
plot(tb,xb(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('population size variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)


%%
%set up vector of control costs for each control- split as control cost and
%cost of model state
figure;
Cost = [costm2,nan,costv2,nan,costb2;costm-costm2,nan,costv-costv2,nan,costb-costb2]';
%plot as stacked bar chart
bar(Cost,'stacked');
hold on;
ax = gca;
ax.XTickLabels = {'Mean','','Variance','','Both'};
grid on;
xlabel('Control Aim')
ylabel('Cost')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
