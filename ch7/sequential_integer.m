% sequentiallly determined integer control on two genotype modelsfor
% sub-intervals of length t_int
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
% set up ode model from functions
dx_m=@(t,x)geno_2_mean_dyn(t,x);

dx_v=@(t,x)geno_2_meta_dyn(t,x);

init_m=[25;5];
init_v=[25;25;10;10];
% run to eqm to get initial conditions
[tm,xm]=ode45(dx_m, [0 15],init_m);
[tv,xv]=ode45(dx_v, [0 150],init_v);
% set cost function parameters
alp1=0.5/(6.23^2);
alpr=0.5/(0.32^2);
MT=20.64;
pT=0.10;
% ste initail conditoins as eqm
init=xm(end,:);
eqm2=xv(end,:);
% set time vector for treatmentchange times
t_int=0.25;

tvec=0:t_int:15;
[~,m]=size(tvec);
% set up control recording vector
uvec=zeros(size(tvec));
cost025=0;
xs=zeros(0,2);
ts=zeros(0,1);
%loop over time vector
for ii=1:m-1
    % MEAN FIELD VERSION - comment out as appropriate
    % run subinterval without treat
    u=0;
    tu=tvec(ii):0.01: tvec(ii+1);
    [~,x1]=ode45(dx_m,tu,init);
    % run subinterval with treatment
    u=1;
    [~,x2]=ode45(dx_m,tu,init);
    % calcualte costs of each treatment
    cost1 = alp1*(x1(:,1)+x1(:,2)-MT).^2...
    ...+alpr*(x1(:,2)./(x1(:,1)+x1(:,2))-pT).^2;
    cost1=0.01*trapz(cost1);
      cost2 = alp1*(x2(:,1)+x2(:,2)-MT).^2...
      ...+alpr*(x2(:,2)./(x2(:,1)+x2(:,2))-pT).^2;
    cost2=0.01*trapz(cost2)+t_int;
    
    % VAriance VERSION - comment out as appropriate
    % run subinterval without treat
%     u=0;
%       init=eqm2;
%     tu=tvec(ii):0.01: tvec(ii+1);
%     [~,x1]=ode45(dx_v,tu,init);
%     % run subinterval with treatment
%     u=1;
%     [~,x2]=ode45(dx_m,tu,init);
%     % calcualte costs of each treatment
%     cost1 = alp1*(x1(:,1)+x1(:,2)-MT).^2...
%...+alpr*(x1(:,2)./(x1(:,1)+x1(:,2))-pT).^2+ alp2*(x1(:,3)+x1(:,4)-VT).^2;
%     cost1=0.01*trapz(cost1);
%       cost2 = alp1*(x2(:,1)+x2(:,2)-MT).^2...
%...+alpr*(x2(:,2)./(x2(:,1)+x2(:,2))-pT).^2+ alp2*(x2(:,3)+x2(:,4)-VT).^2;
%     cost2=0.01*trapz(cost2)+t_int;
    
    ts=horzcat(ts,tu);
    % compare costs and choose best option
    if cost1<cost2
        % record action
        uvec(ii)=0;
        % updat ei.c for next subinterval
        init=x1(end,:);
        % update total cost
        cost025=cost025+cost1;
        xs=vertcat(xs,x1);
    else
        
    uvec(ii)=1;
        init=x2(end,:);
        cost025=cost025+cost2;
           xs=vertcat(xs,x2);
    end
end
uvec(end)=uvec(end-1);
sum025=sum(uvec(1:end-1));
% plot control and response
figure;
subplot(2,1,1)
stairs(tvec,uvec,'k','LineWidth',4)
xlabel('time')
ylabel('control on/off')
subplot(2,1,2)
plot(ts,xs(:,1)+xs(:,2),'k','LineWidth',4)
xlabel('time')
ylabel('M_{total}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

figure;
subplot(2,1,1)
plot(ts,xs(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('M_S')

subplot(2,1,2)
plot(ts,xs(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('M_R')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

%% repeat for different sub-int length
init=xm(end,:);
eqm2=xv(end,:);

t_int=0.5;

tvec=0:t_int:15;
[~,m]=size(tvec);
uvec=zeros(size(tvec));
cost05=0;
xs=zeros(0,2);
ts=zeros(0,1);
for ii=1:m-1
    u=0;
    tu=tvec(ii):0.01: tvec(ii+1);
    [~,x1]=ode45(dx_m,tu,init);
    u=1;
    [~,x2]=ode45(dx_m,tu,init);
    
    cost1 = alp1*(x1(:,1)+x1(:,2)-MT).^2 ...
    ...+alpr*(x1(:,2)./(x1(:,1)+x1(:,2))-pT).^2;
    cost1=0.01*trapz(cost1);
      cost2 = alp1*(x2(:,1)+x2(:,2)-MT).^2 ...
      ...+alpr*(x2(:,2)./(x2(:,1)+x2(:,2))-pT).^2;
    cost2=0.01*trapz(cost2)+t_int;
     % VAriance VERSION - comment out as appropriate
    % run subinterval without treat
%     u=0;
%       init=eqm2;
%     tu=tvec(ii):0.01: tvec(ii+1);
%     [~,x1]=ode45(dx_v,tu,init);
%     % run subinterval with treatment
%     u=1;
%     [~,x2]=ode45(dx_m,tu,init);
%     % calcualte costs of each treatment
%     cost1 = alp1*(x1(:,1)+x1(:,2)-MT).^2...
%...+alpr*(x1(:,2)./(x1(:,1)+x1(:,2))-pT).^2+ alp2*(x1(:,3)+x1(:,4)-VT).^2;
%     cost1=0.01*trapz(cost1);
%       cost2 = alp1*(x2(:,1)+x2(:,2)-MT).^2...
%...+alpr*(x2(:,2)./(x2(:,1)+x2(:,2))-pT).^2+ alp2*(x2(:,3)+x2(:,4)-VT).^2;
%     cost2=0.01*trapz(cost2)+t_int;

    ts=horzcat(ts,tu);
    if cost1<cost2
        uvec(ii)=0;
        init=x1(end,:);
        cost05=cost05+cost1;
        xs=vertcat(xs,x1);
    else
        
    uvec(ii)=1;
        init=x2(end,:);
        cost05=cost05+cost2;
           xs=vertcat(xs,x2);
    end
end
uvec(end)=uvec(end-1);
sum05=sum(uvec(1:end-1));
% plot control nd response
figure;
subplot(2,1,1)
stairs(tvec,uvec,'k','LineWidth',4)
xlabel('time')
ylabel('control on/off')
subplot(2,1,2)
plot(ts,xs(:,1)+xs(:,2),'k','LineWidth',4)
xlabel('time')
ylabel('M_{total}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

figure;
subplot(2,1,1)
plot(ts,xs(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('M_S')

subplot(2,1,2)
plot(ts,xs(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('M_R')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

%% repeat for different sub-int length

init=xm(end,:);
eqm2=xv(end,:);

t_int=1;

tvec=0:t_int:15;
[~,m]=size(tvec);
uvec=zeros(size(tvec));
cost10=0;
xs=zeros(0,2);
ts=zeros(0,1);
for ii=1:m-1
    u=0;
    tu=tvec(ii):0.01: tvec(ii+1);
    [~,x1]=ode45(dx_m,tu,init);
    u=1;
    [~,x2]=ode45(dx_m,tu,init);
    
    cost1 = alp1*(x1(:,1)+x1(:,2)-MT).^2 ...
        ...+alpr*(x1(:,2)./(x1(:,1)+x1(:,2))-pT).^2;
    cost1=0.01*trapz(cost1);
      cost2 = alp1*(x2(:,1)+x2(:,2)-MT).^2 ...
          ...+alpr*(x2(:,2)./(x2(:,1)+x2(:,2))-pT).^2;
    cost2=0.01*trapz(cost2)+t_int;
     % VAriance VERSION - comment out as appropriate
    % run subinterval without treat
%     u=0;
%       init=eqm2;
%     tu=tvec(ii):0.01: tvec(ii+1);
%     [~,x1]=ode45(dx_v,tu,init);
%     % run subinterval with treatment
%     u=1;
%     [~,x2]=ode45(dx_m,tu,init);
%     % calcualte costs of each treatment
%     cost1 = alp1*(x1(:,1)+x1(:,2)-MT).^2 ...
% ...+alpr*(x1(:,2)./(x1(:,1)+x1(:,2))-pT).^2+ alp2*(x1(:,3)+x1(:,4)-VT).^2;
%     cost1=0.01*trapz(cost1);
%       cost2 = alp1*(x2(:,1)+x2(:,2)-MT).^2 ...
% ...+alpr*(x2(:,2)./(x2(:,1)+x2(:,2))-pT).^2+ alp2*(x2(:,3)+x2(:,4)-VT).^2;
%     cost2=0.01*trapz(cost2)+t_int;

    ts=horzcat(ts,tu);
    if cost1<cost2
        uvec(ii)=0;
        init=x1(end,:);
        cost10=cost10+cost1;
        xs=vertcat(xs,x1);
    else
        
    uvec(ii)=1;
        init=x2(end,:);
        cost10=cost10+cost2;
           xs=vertcat(xs,x2);
    end
end
uvec(end)=uvec(end-1);
sum10=sum(uvec(1:end-1));
figure;
subplot(2,1,1)
stairs(tvec,uvec,'k','LineWidth',4)
xlabel('time')
ylabel('control on/off')
subplot(2,1,2)
plot(ts,xs(:,1)+xs(:,2),'k','LineWidth',4)
xlabel('time')
ylabel('M_{total}')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

figure;
subplot(2,1,1)
plot(ts,xs(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('M_S')

subplot(2,1,2)
plot(ts,xs(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('M_R')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
% set up cost vector 
cost_v = [cost025-0.25*sum025,sum025*0.25;nan,nan;...
    ...cost05-0.5*sum05,sum05*0.5;nan,nan;cost10-1*sum10,sum10*1];
% plot cost comparison
figure;

bar(cost_v,'stacked');
hold on;
ax = gca;
ax.XTickLabels = {'\tau=0.25','','\tau=0.5','','\tau=1'};
grid on;
xlabel('Discrete Interval length')
ylabel('Cost')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
