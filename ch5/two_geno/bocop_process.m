% code to take results of bocop optimisation and plot them for two genotype
% model
% set parameter values
clear all
close all
global beta mu_M D_M psi1  u1 tu
beta=1.5;
mu_M= 1/(5*365);
D_M=1/20;

psi1=0.9;
% set up ode model with defined functions
dx_m=@(t,x)geno_2_mean_dyn_t(t,x);
dx_v=@(t,x)geno_2_meta_dyn_t(t,x);
% ste initial condition = to bocop i.c's
init_m = [29.6;3.5];
init_v=[29.2;2.8;14.9;2.6];

%import bocop data
bocop1 = readtable('thwo_geno_mods.xlsx');

% extract controls basedon v,m,v&m
 u_mm=bocop1{1:end,{'param2mean'}};
 u_vm=bocop1{1:end,{'param2var'}};
  u_vv=bocop1{1:end,{'param4var'}};
   u_vb=bocop1{1:end,{'param6var'}};
% u_v=bocop1{1:end,{'var'}};
% 
% u_mm=bocop1{1:end,{'meanmod'}};
% 
% set up time vector for controls
m=length(u_vm);
tu=0:30/(m-1):30;

% set control as each different controland run on corresponding simulation
u1=u_mm;

[tm,xm]=ode45(dx_m,[0 30],init_m);

u1=u_vm;
[tvm,xvm]=ode45(dx_v,[0 30],init_v);

u1=u_vv;
[tvv,xvv]=ode45(dx_v,[0 30],init_v);

u1=u_vb;
[tvb,xvb]=ode45(dx_v,[0 30],init_v);

% plotresulst
figure;
subplot(2,1,1)
plot(tu,u_mm,'k','LineWidth',4)
xlabel('Time')
ylabel('control intensity')
subplot(2,1,2)
plot(tm,xm(:,1)+xm(:,2),'k','LineWidth',4)
xlabel('Time')
ylabel('M_{total}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
figure;
subplot(2,1,1)
plot(tm,xm(:,1),'r','LineWidth',4)
xlabel('Time')
ylabel('M_{s}')
subplot(2,1,2)
plot(tm,xm(:,2),'b','LineWidth',4)
xlabel('Time')
ylabel('M_{r}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

figure;
subplot(3,1,1)
plot(tu,u_vm,'k','LineWidth',4)
xlabel('Time')
ylabel('control intensity')
subplot(3,1,2)
plot(tvm,xvm(:,1)+xvm(:,2),'r','LineWidth',4)
xlabel('Time')
ylabel('M_{total}')
subplot(3,1,3)
plot(tvm,xvm(:,3)+xvm(:,4),'b','LineWidth',4)
xlabel('Time')
ylabel('V_{total}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
figure;
subplot(2,2,1)
plot(tvm,xvm(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('susceptible mean')
subplot(2,2,2)
plot(tvm,xvm(:,2),'--r','LineWidth',4)
xlabel('time')
ylabel('resistant mean')
subplot(2,2,3)
plot(tvm,xvm(:,3),'b','LineWidth',4)
xlabel('time')
ylabel('susceptible variance')
subplot(2,2,4)
plot(tvm,xvm(:,4),'--b','LineWidth',4)
xlabel('time')
ylabel('resistant variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)


figure;
subplot(3,1,1)
plot(tu,u_vv,'k','LineWidth',4)
xlabel('Time')
ylabel('control intensity')
subplot(3,1,2)
plot(tvv,xvv(:,1)+xvv(:,2),'r','LineWidth',4)
xlabel('Time')
ylabel('M_{total}')
subplot(3,1,3)
plot(tvv,xvv(:,3)+xvv(:,4),'b','LineWidth',4)
xlabel('Time')
ylabel('V_{total}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
figure;
subplot(2,2,1)
plot(tvv,xvv(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('susceptible mean')
subplot(2,2,2)
plot(tvv,xvv(:,2),'--r','LineWidth',4)
xlabel('time')
ylabel('resistant mean')
subplot(2,2,3)
plot(tvv,xvv(:,3),'b','LineWidth',4)
xlabel('time')
ylabel('susceptible variance')
subplot(2,2,4)
plot(tvv,xvv(:,4),'--b','LineWidth',4)
xlabel('time')
ylabel('resistant variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)


figure;
subplot(3,1,1)
plot(tu,u_vb,'k','LineWidth',4)
xlabel('Time')
ylabel('control intensity')
subplot(3,1,2)
plot(tvb,xvb(:,1)+xvb(:,2),'r','LineWidth',4)
xlabel('Time')
ylabel('M_{total}')
subplot(3,1,3)
plot(tvb,xvb(:,3)+xvb(:,4),'b','LineWidth',4)
xlabel('Time')
ylabel('V_{total}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
figure;
subplot(2,2,1)
plot(tvb,xvb(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('susceptible mean')
subplot(2,2,2)
plot(tvb,xvb(:,2),'--r','LineWidth',4)
xlabel('time')
ylabel('resistant mean')
subplot(2,2,3)
plot(tvb,xvb(:,3),'b','LineWidth',4)
xlabel('time')
ylabel('susceptible variance')
subplot(2,2,4)
plot(tvb,xvb(:,4),'--b','LineWidth',4)
xlabel('time')
ylabel('resistant variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

