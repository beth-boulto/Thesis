%% code to import data from the bocop optimisation of control for the mean 
% and variance models and plot the figures
clear all
close all
%set parameters
%% bocop plotting code for mean and metapop models
% set parameter
global beta mu_M D_M tu u
beta=1.5;
mu_M=1/(5*365);
D_M=1/20;
u=0;
% set up ode models using defined functions from previous
dx=@(t,x)metapop(t,x);
dxm=@(t,x)mean_field(t,x);
%import bocop data
bocop1 = readtable('corrected_basic_model.xlsx');

% extract controls basedon v,m,v&m
u_m=bocop1{1:end,{'mean'}};
u_vm=bocop1{1:end,{'meanvar'}};
u_v=bocop1{1:end,{'var'}};

u_mm=bocop1{1:end,{'meanmod'}};

m=length(u_vm);
tu=0:10/(m-1):10;

%% mean control
u=u_m;
% run control with mean control imported from bocop
[tm,xm]=ode45(dx,[0 10],[29.5,15]);
% calculate cost from cost function
Y = 0.00908*(xm(:,1)-9).^2;
costm= trapz(tu,u.^2);
costm=costm+trapz(tm,Y);
costm2 =  trapz(tu,u.^2);
% plot control
figure;
subplot(3,1,1)
plot(tu,u,'k','LineWidth',4)
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


%% variance control - only variance is penalised
% set control as imported from bocop
u=u_v;
% simulate on metapop model
[tv,xv]=ode45(dx,[0 10],[29.5,15]);
% calculate cost via cost function
Y = 0.0626*(xv(:,2)-7).^2;
costv= trapz(tu,u.^2);
costv=costv+trapz(tv,Y);
costv2 =  trapz(tu,u.^3);

% plot control and response
figure;
subplot(3,1,1)
plot(tu,u,'k','LineWidth',4)
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
% set control as bocop import
u=u_vm;
% simulate ressposnse to control
[tb,xb]=ode45(dx,[0 10],[29.5,15]);
% calculate cost via cost function
Y = 0.00454*(xb(:,1)-9).^2+0.0313*(xb(:,2)-7).^2;
costb= trapz(tu,u.^2);
costb=costb+trapz(tb,Y);
costb2 =  trapz(tu,u.^2);
% plot control and response
figure;
subplot(3,1,1)
plot(tu,u,'k','LineWidth',4)
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

%% mean model mean control
%set control to bocop result
u=u_mm;
% simulate response to control
[tmm,xmm]=ode45(dxm,[0 10],30);
% calculate cost of control
Y = 1/36*(xmm-9.5).^2;
costmm=trapz(tu,u);
costmm=costmm+trapz(tmm,Y);
% plot control and response
figure;
subplot(2,1,1)
plot(tu,u,'k','LineWidth',4)
xlabel('time')
ylabel('control intensity')
subplot(2,1,2)
plot(tmm,xmm(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('mean population size')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
%% cost comparison
figure;
% set up vector of all costs of control - as control cost, cost due to state
Cost = [costm2,nan,costv2,nan,costb2;costm-costm2,nan,costv-costv2,nan,costb-costb2]';

% plot as stacked bar chart
bar(Cost,'stacked');
hold on;
ax = gca;
ax.XTickLabels = {'Mean','','Variance','','Both'};
grid on;
xlabel('Control Aim')
ylabel('Cost')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
