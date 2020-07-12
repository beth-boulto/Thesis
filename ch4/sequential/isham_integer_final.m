%% code to calculate control on the ISham model on sub-intervals sequentially
% with sub-intervals of length tau1. The code then plots the results and
% calculates the total cost of the control.
clear all
close all
global phi alpha mu_M p r u1 tu
% set parameter values
phi = 52/12;
mu_H = 0;
mu_M1= 1/12;
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
  dx= @(t,x)isham_treat(t,x);

tau1=1;
% set treatment time vector
tu=0:tau1:10;

m=length(tu);
ua=zeros(1,m);

% set imitial state from untreatedd eqm
init1=[25,45];
xt=zeros(2,0);
ts=zeros(0,1);
% set cost function parameters
alp1=0/(2*11.51^2);
alp2=2/(2*20.755^2);
MT=1.98;
VT=3.49;
% set initial cost
costt=0;
for ii=1:m-1
    % set time vector for curent treatment interval
   tt=tu(ii):0.01:tu(ii+1);
   % simulate interval without treatment
   mu_M=mu_M1;
   [~,x1]=ode45(dx,tt, init1);
    % simulate interval with treatment
   mu_M=mu_M1+1;
   [~,x2]=ode45(dx,tt,init1);
   % calculate cost functions for both scenarios
   Y1=alp1*(x1(:,1)-MT).^2+alp2*(x1(:,2)-VT).^2;
   Y1=0.01*trapz(Y1);
   cost1=Y1;
   Y1=alp1*(x2(:,1)-MT).^2+alp2*(x2(:,2)-VT).^2;
   Y1=0.01*trapz(Y1);
   cost2=Y1+1*(tu(ii+1)-tu(ii));
   ts=horzcat(ts,tt);
   % compare costs - choose lowest
   if cost1<cost2
       % set initial conditions of next interval using end conditions of
       % simulation corresponding to chosen action
       ua(ii)=0;
       init1=x1(end,:);
       xt=vertcat(xt,x1);
       % update total cost
       costt=costt+cost1;
   else
       ua(ii)=1;
       init1=x2(end,:);
            xt=vertcat(xt,x2);
            costt=costt+cost2;
   end
   
    
end
% plot final control and response 
u=ua;
figure;
subplot(3,1,1)
stairs(tu,u,'k','LineWidth',4)
xlabel('time')
ylabel('control on/off')
subplot(3,1,2)
plot(ts,xt(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('mean population size')
subplot(3,1,3)
plot(ts,xt(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('population size variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
%% repeat with differnt treatment interval length
tu=0:0.5:10;

m=length(tu);
u2=zeros(1,m);
init1=[25,45];
xt=zeros(2,0);
ts=zeros(0,1);

costt2=0;
for ii=1:m-1
   tt=tu(ii):0.01:tu(ii+1);
   mu_M=mu_M1;
   [~,x1]=ode45(dx,tt, init1);
   mu_M=mu_M1+1;
   [~,x2]=ode45(dx,tt,init1);
   
   Y1=alp1*(x1(:,1)-MT).^2+alp2*(x1(:,2)-VT).^2;
   Y1=0.01*trapz(Y1);
   cost1=Y1;
   Y1=alp1*(x2(:,1)-MT).^2+alp2*(x2(:,2)-VT).^2;
   Y1=0.01*trapz(Y1);
   cost2=Y1+1*(tu(ii+1)-tu(ii));
   ts=horzcat(ts,tt);
   if cost1<cost2
       u2(ii)=0;
       init1=x1(end,:);
       xt=vertcat(xt,x1);
       costt2=costt2+cost1;
   else
       u2(ii)=1;
       init1=x2(end,:);
            xt=vertcat(xt,x2);
            costt2=costt2+cost2;
   end
   
    
end
u2(end)=u2(end-1);
u=u2;
figure;
subplot(3,1,1)
stairs(tu,u,'k','LineWidth',4)
xlabel('time')
ylabel('control on/off')
subplot(3,1,2)
plot(ts,xt(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('mean population size')
subplot(3,1,3)
plot(ts,xt(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('population size variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
%% repeat with differnt treatment interval length

tu=0:0.25:10;

m=length(tu);
u3=zeros(1,m);
init1=[25,45];
xt=zeros(2,0);
ts=zeros(0,1);


costt3=0;
for ii=1:m-1
   tt=tu(ii):0.01:tu(ii+1);
   mu_M=mu_M1;
   [~,x1]=ode45(dx,tt, init1);
   mu_M=mu_M1+1;
   [~,x2]=ode45(dx,tt,init1);
   
   Y1=alp1*(x1(:,1)-MT).^2+alp2*(x1(:,2)-VT).^2;
   Y1=0.01*trapz(Y1);
   cost1=Y1;
   Y1=alp1*(x2(:,1)-MT).^2+alp2*(x2(:,2)-VT).^2;
   Y1=0.01*trapz(Y1);
   cost2=Y1+1*(tu(ii+1)-tu(ii));
   ts=horzcat(ts,tt);
   if cost1<cost2
       u3(ii)=0;
       init1=x1(end,:);
       xt=vertcat(xt,x1);
       costt3=costt3+cost1;
   else
       u3(ii)=1;
       init1=x2(end,:);
            xt=vertcat(xt,x2);
            costt3=costt3+cost2;
   end
   
    
end

u3(end)=u3(end-1);
u=u3;
figure;
subplot(3,1,1)
stairs(tu,u,'k','LineWidth',4)
xlabel('time')
ylabel('control on/off')
subplot(3,1,2)
plot(ts,xt(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('mean population size')
subplot(3,1,3)
plot(ts,xt(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('population size variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
%% repeat with differnt treatment interval length
 
tu=0:5:10;

m=length(tu);
u4=zeros(1,m);
init1=[25,45];
xt=zeros(2,0);
ts=zeros(0,1);


costt4=0;
for ii=1:m-1
   tt=tu(ii):0.01:tu(ii+1);
   mu_M=mu_M1;
   [~,x1]=ode45(dx,tt, init1);
   mu_M=mu_M1+1;
   [~,x2]=ode45(dx,tt,init1);
   
   Y1=alp1*(x1(:,1)-MT).^2+alp2*(x1(:,2)-VT).^2;
   Y1=0.01*trapz(Y1);
   cost1=Y1;
   Y1=alp1*(x2(:,1)-MT).^2+alp2*(x2(:,2)-VT).^2;
   Y1=0.01*trapz(Y1);
   cost2=Y1+1*(tu(ii+1)-tu(ii));
   ts=horzcat(ts,tt);
   if cost1<cost2
       u4(ii)=0;
       init1=x1(end,:);
       xt=vertcat(xt,x1);
       costt4=costt4+cost1;
   else
       u4(ii)=1;
       init1=x2(end,:);
            xt=vertcat(xt,x2);
            costt4=costt+cost2;
   end
   
    
end

u4(end)=u4(end-1);
u=u4;
figure;
subplot(3,1,1)
stairs(tu,u,'k','LineWidth',4)
xlabel('time')
ylabel('control on/off')
subplot(3,1,2)
plot(ts,xt(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('mean population size')
subplot(3,1,3)
plot(ts,xt(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('population size variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
%% run optimal control for corresponding cost function using bocp import
bocop1 = readtable('isham.xlsx');
u_m = bocop1{1:end,{'mean_11_51var_20_755'}};
u_vm=bocop1{1:end,{'meanvar'}};
u_v=bocop1{1:end,{'Var2'}};
m=length(u_vm);
tu=0:10/(m-1):10;
dx= @(t,x)isham_treat(t,x);

u1=u_m;
[tm,xm]=ode45(dx,[0 10],[25,45]);
% calculate cost of optimal control
Y = alp1*(xm(:,1)-1.98).^2+alp2*(xm(:,2)-VT).^2;
costm=trapz(tm,Y);
costm2=trapz(tu,u1);

% set up cost vectors to compare costs- plot as bar graph
figure;
cost_vec=[costt4-sum(u4(1:end-1))*5,sum(u4(1:end-1))*5;nan,nan;costt-sum(ua(1:end-1))*1,sum(ua(1:end-1))*1;nan,nan;costt2-sum(u2(1:end-1))*0.5,sum(u2(1:end-1))*0.5;nan,nan;costt3-0.25*sum(u3(1:end-1)),sum(u3(1:end-1))*0.25;nan,nan;costm, costm2];
bar(cost_vec,'stacked');
hold on;
ax = gca;
ax.XTickLabels = {'\tau=5','','\tau=1','','\tau=0.5','','\tau=0.25','','Continuous'};
grid on;
xlabel('Discrete Interval length')
ylabel('Cost')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

