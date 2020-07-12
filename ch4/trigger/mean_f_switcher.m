%% code to simulate switching of treatment on mean field model when 
%specific conditions are met using a switch system. The codes then plots
%the results and calculates the cost of control.
clear all 
close all
% set parameter values
global beta mu_M D_M MT T_start
beta=1.5;
D_M=1/20;
mu_M1=1/(5*365);
% set trigger conditions
alp1=1/(11.74^2);
MT=9.5+11.74;
% set initial condition
init1=[30];
% set ode model using defined funtcion
dx=@(t,x)mean_field(t,x);
% set final time
T_f=10;
% set current/initial  time
t_s=0;
% set recording vector for treatment times
tc=zeros(1,0);
uc=zeros(0,1);
% set initial condiions
init=init1;
ts=zeros(0,1);
xs=zeros(0,1);
% while current time less than final time
while t_s<T_f
         tc=horzcat(tc,t_s);
         % fix interval start at current time
         T_start=t_s;
   uc=horzcat(uc,0);
   % begin untreated period
    mu_M=mu_M1;
    % set event function to stop simulation under set conditions
    optfun1=@(t,x)option_fun_mean_on(t,x);
    opts=odeset('Events',optfun1);
    % run simulation
    [t1,x1]=ode45(dx,[t_s T_f],init,opts);
    
   ts=vertcat(ts,t1);
   xs=vertcat(xs,x1);
   % update current time
   t_s=ts(end);
   init=xs(end,:);
   % if ended before final end time 
   if t_s<T_f
       % set interval start time as current time
%        T_start=t_s;
        tc=horzcat(tc,t_s);
    uc=horzcat(uc,1);
    % set parasite death for treatment
        mu_M=mu_M1+1;
        % set optionfunction to stop treatment under chosen conditions
%         optfun2=@(t,x)option_fun_mean_off(t,x);
%     opts=odeset('Events',optfun2);
%     [t1,x1]=ode45(dx,[t_s T_f],init,opts);
% if treatment is set on timer
t_e=min(T_f, t_s+0.4);
% run simulation
[t1,x1]=ode45(dx,[t_s t_e],init);

    ts=vertcat(ts,t1);
   xs=vertcat(xs,x1);
   %update current time and initial conditions of next simulation
   t_s=ts(end);
   init=xs(end,:);
   end
   
end
tc=horzcat(tc,10);
uc=horzcat(uc,0);
% plotcontrol and resposne 
figure;
subplot(2,1,1)
stairs(tc,uc,'k','LineWidth',4)
xlabel('time')
ylabel('control on/off')
subplot(2,1,2)
plot(ts,xs,'r','LineWidth',4)
xlabel('time')
ylabel('mean')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
% calculate cost using cost function from previous
MT=9.5;
cr1=alp1*(xs-MT).^2;
cr1=trapz(ts,cr1);
m=length(uc);
cu1=0;
for ii=1:m-1
    cu1=cu1+uc(ii)*(tc(ii+1)-tc(ii));
end
% bocop comparison to continuous optimal control
bocop_b=readtable('bocop_cost.xlsx');
crb=bocop_b{1:end,{'meanmodr'}};
cru=bocop_b{1:end,{'meanmodu'}};
% plot cost comparison
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

   