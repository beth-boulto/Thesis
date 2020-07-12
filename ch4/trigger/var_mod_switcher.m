%% variance model switching system code which simulates the variance 
%inclusive model with and turns the control on or off when the switch is
%activtaed. The code also plots the resulst of the treatment and calculates
%the cost of control.
clear all 
close all
% set parameters
global beta mu_M D_M MT VT T_start
beta=1.5;
D_M=1/20;
mu_M1=1/(5*365);
% set cost function parameters
alp1=0.5/(10.25^2);
alp2=0.5/4^2;
% set condition values for triggers
MT=9+10.25;
VT=7+4;
% set initial condition as approx untreated eqm
init1=[29.5;15];
% set ode model
dx=@(t,x)metapop(t,x);
% set final and initial time
T_f=10;
t_s=0;
tc=zeros(1,0);
uc=zeros(0,1);
init=init1;
ts=zeros(0,1);
xs=zeros(0,1);
% while current time less than final time
while t_s<T_f
    % set start of untrreated period
         tc=horzcat(tc,t_s);
         T_start=t_s;
   uc=horzcat(uc,0);
   % set parasite death rate ( no treat)
    mu_M=mu_M1;
    % choose switching function
    
    %optfun1=@(t,x)option_fun_var_var_on(t,x);
     %optfun1=@(t,x)option_fun_var_var_on(t,x);
      optfun1=@(t,x)option_fun_var_mean_on(t,x);
      
      % set ode options to stop when evet occurs according to switch
      % function
    opts=odeset('Events',optfun1);
    % simulate until either final time or event occurs
    [t1,x1]=ode45(dx,[t_s T_f],init,opts);
   ts=vertcat(ts,t1);
   xs=vertcat(xs,x1);
   t_s=ts(end);
   init=xs(end,:);
   
   % if current time less than fimal time 
   if t_s<T_f
       % set start time of treated section
       T_start=t_s;
       tc=horzcat(tc,t_s);
   uc=horzcat(uc,1);
   % set parasite death plus treatment
       mu_M=mu_M1+1;
       
       % set option function for stopping 
        %optfun2=@(t,x)option_fun_var_var_off(t,x);
        %optfun2=@(t,x)option_fun_var_var_off(t,x);
%         optfun2=@(t,x)option_fun_var_both_off(t,x);
% set ode options to stop when event occurs
%     opts=odeset('Events',optfun2);
%     [t1,x1]=ode45(dx,[t_s T_f],init,opts);
% if choesen set end of treatment time
t_e=min(T_f,t_s+0.4);
% simulate until treatment stops
[t1,x1]=ode45(dx,[t_s t_e],init);
% concatenate with previosu intervals
    ts=vertcat(ts,t1);
   xs=vertcat(xs,x1);
   % update current ime and initial conditions
   t_s=ts(end);
   init=xs(end,:);
   end
   
end
tc=horzcat(tc,10);
uc=horzcat(uc,0);
% plot control and response 
figure;
subplot(3,1,1)
stairs(tc,uc,'k','LineWidth',4)
xlabel('time')
ylabel('control on/off')
subplot(3,1,2)
plot(ts,xs(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('mean')
subplot(3,1,3)
plot(ts,xs(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
% calculate cost of control
MT=9;
VT=7;
cr1=alp1*(xs(:,1)-MT).^2+alp2*(xs(:,2)-VT).^2;
cr1=trapz(ts,cr1);
cu1=0;
m=length(uc);
for ii=1:m-1
    cu1=cu1+uc(ii)*(tc(ii+1)-tc(ii));
end

% compare with cts optimal control
bocop_b=readtable('bocop_cost.xlsx');
crb=bocop_b{1:end,{'vmodbr'}};
cru=bocop_b{1:end,{'vmodbu'}};
% crb=bocop_b{1:end,{'vmodvr'}};
% cru=bocop_b{1:end,{'vmodvu'}};
% crb=bocop_b{1:end,{'vmodbr'}};
% cru=bocop_b{1:end,{'vmodbu'}};
figure;
cost_v=[cr1,cu1;nan,nan;crb,cru];
bar(cost_v,'stacked');
hold on;
ax = gca;
ax.XTickLabels = {'Switched','','Continuous optimal'};
grid on;
xlabel('Method')
ylabel('Cost')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

   