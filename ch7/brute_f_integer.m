%% code to do brute force optimisation on 2 genotype var inclusive model 
%or mean field model with sub-interval length tau1

clear all
close all
% set parameters
global beta mu_M D_M psi1  u1 tu
beta=1.5;
mu_M= 1/(5*365);
D_M=1/20;
u=0;
psi1=0.9;
% set treatment at zero
tu=[0 150];
u1=[0,0];
% set ode models with predefined fucntions
dx_m=@(t,x)geno_2_mean_dyn_t(t,x);
dx_v=@(t,x)geno_2_meta_dyn_t(t,x);
init_m=[25;5];
init_v=[25;25;10;10];
% run to eqm
[tm,xm]=ode45(dx_m, [0 150],init_m);
[tv,xv]=ode45(dx_v, [0 150],init_v);
% set inital conditions as eqm
init_m=xm(end,:);
init_v=xv(end,:);
% set start and end times
t_0=0;
T_f=20;
% set up vector of treatment change times
tau1=T_f/10;
tu=t_0:tau1:T_f;
m=length(tu);
% find all combos of treat no treat for time vector
u_mat=de2bi(0:2^(m)-1);
[n,~]=size(u_mat);
t_l=tu(2)-tu(1);
xt=zeros(4,0);
ts=zeros(0,1);
%cost3=zeros(1,2^m);
% set initial cost to inf for comparison
cost=inf;
treat_num=0;
% set cost function parameters
alp1=0;
alp2=1/((1.73)^2);
MT= 19;
VT = 14;
% loop over all treatment combos
for ii=1:n
    
     u1=u_mat(ii,:);
     %VARIANCE MODEL VERSION - choose intended model
     % run for chosen treatment
     [ts,xs]=ode45(dx_v,[0 T_f], init_v);
     % calculate cost
    Y1=alp1*(xs(:,1)+xs(:,2)-MT).^2+alp2*(xs(:,3)+xs(:,4)-VT).^2;
 % MEAN FIELD MODEL VERSION   - choose intended model
%     [ts,xs]=ode45(dx_m,[0 T_f], init_m);
%      calculate cost
%     Y1=alp1*(xs(:,1)+xs(:,2)-MT).^2;
    
    cost2 = sum(u)*t_l+trapz(ts,Y1);
    
    %cost2=cost2-u(end)*1;
    % compare cost with current best
    if cost2<cost
        % if better update current cost and treatment choice
        cost=cost2;
        treat_num=ii;
        t1=ts;
        x1=xs;
    end
    
    
end
% plot chosen treatment and response - variance model, adapt for mean model
% by removing x1(:,3) and x1(:,4) options
figure;
subplot(3,1,1)
stairs(tu,u_mat(treat_num,:),'k','LineWidth',4)
xlabel('time')
ylabel('treatment on/off')
subplot(3,1,2)
plot(t1,x1(:,1)+x1(:,2),'r','LineWidth',4)
xlabel('time')
ylabel('M_{total}')
subplot(3,1,3)
plot(t1,x1(:,3)+x1(:,4),'b','LineWidth',4)
xlabel('time')
ylabel('V_{total}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)


figure;
subplot(2,2,1)
plot(t1,x1(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('M_{s}')
subplot(2,2,2)
plot(t1,x1(:,2),'--r','LineWidth',4)
xlabel('time')
ylabel('M_{r}')
subplot(2,2,3)
plot(t1,x1(:,3),'b','LineWidth',4)
xlabel('time')
ylabel('V_{s}')
subplot(2,2,4)
plot(t1,x1(:,4),'--b','LineWidth',4)
xlabel('time')
ylabel('V_{r}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

