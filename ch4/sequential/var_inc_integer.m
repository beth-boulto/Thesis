%% code to calculate control on the mean field and variance inclusive model
%on sub-intervals sequentially
% with sub-intervals of length tau1. The code then plots the results and
% calculates the total cost of the control.
%set parameters
global beta mu_M D_M tu u
beta=1.5;
mu_M1=1/(5*365);
D_M=1/20;
tsum=0:0.01:10;
% set ode models from functions
dx=@(t,x)metapop(t,x);
dxm=@(t,x)mean_field(t,x);
dx_2=@(t,x)alt_metapop(t,x);
% set when treatment is swicthed on or off
tau1=1;
tu=0:tau1:10;

m=length(tu);
u1=zeros(1,m);
% set initial conditions according to untreatde eqm 
init1=[29.5;15];
xt=zeros(2,0);
ts=zeros(0,1);
% set cost function parameters
alp1=0/(10.25^2);
alp2=1/16;
MT=9;
VT=7;
% set initial cost
cost_t=0;
for ii=1:m-1
    % set time vector for treatment interval
   tt=tu(ii):0.01:tu(ii+1);
   % simulate without treatmemnt
   mu_M=mu_M1;
   [~,x1]=ode45(dx,tt, init1);
   % simulate with treatment
   mu_M=mu_M1+1;
   [~,x2]=ode45(dx,tt,init1);
   % calculate cost without treatment
   Y1=alp1*(x1(:,1)-MT).^2+alp2*(x1(:,2)-VT).^2;
   Y1=trapz(tt,Y1);
   cost1=Y1;
   % calculate cost with treatment applied
   Y1=alp1*(x2(:,1)-MT).^2+alp2*(x2(:,2)-VT).^2;
   Y1=trapz(tt,Y1);
   cost2=Y1+1*(tu(ii+1)-tu(ii));
   ts=horzcat(ts,tt);
   % choose cheapest option
   if cost1<cost2 % no treatment cheaper
       u1(ii)=0;
       % set initial condition of next interval as untreated simulation
       % state at end of interval
       init1=x1(end,:);
       xt=vertcat(xt,x1);
       % update total cost
       cost_t=cost_t+cost1;
   else % treatment cheaper
       u1(ii)=1;
       % set initial condition as state of model at end of treatment 
       init1=x2(end,:);
            xt=vertcat(xt,x2);
            % update total cost
            cost_t=cost_t+cost2;
   end
   
    
end
% finish treatment vector
u1(end)=u1(end-1);
u=u1;
 % plot treatment vector and resposnse to treatment  
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

%% repeat for alternative treatment interval length
tu=0:0.25:10;

m=length(tu);
u3=zeros(1,m);
init1=[29.5;15];
xt=zeros(2,0);
ts=zeros(0,1);

MT=9;
VT=7;
cost_t3=0;
for ii=1:m-1
   tt=tu(ii):0.01:tu(ii+1);
   mu_M=mu_M1;
   [~,x1]=ode45(dx,tt, init1);
   mu_M=mu_M1+1;
   [~,x2]=ode45(dx,tt,init1);
   
   Y1=alp1*(x1(:,1)-MT).^2+alp2*(x1(:,2)-VT).^2;
   Y1=trapz(tt,Y1);
   cost1=Y1;
   Y1=alp1*(x2(:,1)-MT).^2+alp2*(x2(:,2)-VT).^2;
   Y1=trapz(tt,Y1);
   cost2=Y1+1*(tu(ii+1)-tu(ii));
   ts=horzcat(ts,tt);
   if cost1<cost2
       u3(ii)=0;
       init1=x1(end,:);
       xt=vertcat(xt,x1);
       cost_t3=cost_t3+cost1;
   else
       u3(ii)=1;
       init1=x2(end,:);
            xt=vertcat(xt,x2);
            cost_t3=cost_t3+cost2;
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

%% repeat for another treatment length
tu=0:0.5:10;

m=length(tu);
u2=zeros(1,m);
init1=[29.5;15];
xt=zeros(2,0);
ts=zeros(0,1);

MT=9;
VT=7;
cost_t2=0;
for ii=1:m-1
   tt=tu(ii):0.01:tu(ii+1);
   mu_M=mu_M1;
   [~,x1]=ode45(dx,tt, init1);
   mu_M=mu_M1+1;
   [~,x2]=ode45(dx,tt,init1);
   
   Y1=alp1*(x1(:,1)-MT).^2+alp2*(x1(:,2)-VT).^2;
   Y1=trapz(tt,Y1);
   cost1=Y1;
   Y1=alp1*(x2(:,1)-MT).^2+alp2*(x2(:,2)-VT).^2;
   Y1=trapz(tt,Y1);
   cost2=Y1+1*(tu(ii+1)-tu(ii));
   ts=horzcat(ts,tt);
   if cost1<cost2
       u2(ii)=0;
       init1=x1(end,:);
       xt=vertcat(xt,x1);
       cost_t2=cost_t2+cost1;
   else
       u2(ii)=1;
       init1=x2(end,:);
            xt=vertcat(xt,x2);
            cost_t2=cost_t2+cost2;
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

%% repeat for final interval length
tu=0:5:10;

m=length(tu);
u4=zeros(1,m);
init1=[29.5;15];
xt=zeros(2,0);
ts=zeros(0,1);

MT=9;
VT=7;
cost_t4=0;
for ii=1:m-1
   tt=tu(ii):0.01:tu(ii+1);
   mu_M=mu_M1;
   [~,x1]=ode45(dx,tt, init1);
   mu_M=mu_M1+1;
   [~,x2]=ode45(dx,tt,init1);
   
   Y1=alp1*(x1(:,1)-MT).^2+alp2*(x1(:,2)-VT).^2;
   Y1=trapz(tt,Y1);
   cost1=Y1;
   Y1=alp1*(x2(:,1)-MT).^2+alp2*(x2(:,2)-VT).^2;
   Y1=trapz(tt,Y1);
   cost2=Y1+1*(tu(ii+1)-tu(ii));
   ts=horzcat(ts,tt);
   if cost1<cost2
       u4(ii)=0;
       init1=x1(end,:);
       xt=vertcat(xt,x1);
       cost_t4=cost_t4+cost1;
   else
       u4(ii)=1;
       init1=x2(end,:);
            xt=vertcat(xt,x2);
            cost_t4=cost_t4+cost2;
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

%% %import bocop data
bocop1 = readtable('var_mod.xlsx');

% extract controls basedon v,m,v&m
u_m=bocop1{1:end,{'mean_10_25'}};
u_vm=bocop1{1:end,{'mean_var_10_25_4'}};
u_v=bocop1{1:end,{'var_4'}};
% set control to control corresponding to cost function used
u=u_m;
m=length(u_m);
tu=0:10/(m-1):10;
% simulate over same overall time
[tm,xm]=ode45(dx,[0 10],[29.5,15]);
% calculate cost of optimal control
Y =alp1*(xm(:,1)-9).^2+alp2*(xm(:,2)-VT).^2;
costm= trapz(tu,u.^2);
costm2=trapz(tm,Y);
% set up cost vector of integer controls and optimal for same cost function
figure;
cost_vec=[cost_t4-sum(u4(1:end-1))*5,sum(u4(1:end-1))*5;nan,nan; cost_t-sum(u1(1:end-1))*1,sum(u1(1:end-1))*1;nan,nan;cost_t2-sum(u2(1:end-1))*0.5,sum(u2(1:end-1))*0.5;nan,nan;cost_t3-0.25*sum(u3(1:end-1)),sum(u3(1:end-1))*0.25;nan,nan;costm2, costm];
% plot costs as bar chart
bar(cost_vec,'stacked');
hold on;
ax = gca;
ax.XTickLabels = {'\tau=5','','\tau=1','','\tau=0.5','','\tau=0.25','','Continuous'};
grid on;
xlabel('Discrete Interval length')
ylabel('Cost')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
legend('state','application of control');

