%% Code to simulate Isham model with a treatment on most infected hosts for
% a set time using gillespie algorithm. The code then performs a simulation
% of the untreated model using the deterministic moment closure model and
% then calculate the mean and variance of the treated/untreated proportion
% of the population from truncated distributions which are then used to
% simualte the mean/variance of the two separate hosts populations. This
% results are then combined to find the total mean/variance and compare it
% to the stochastic simulation.
close all
clear all
global r p
psi = 52;
mu_H = 0;
mu_M1= 10;
alpha=0.02;
M = 10;
u=20;
    r=0.5;    
C_num=50;

    V=M.^2/r-M;
    
  prob=M/V;
  p=1-prob;

% set up clump distribution
pd1 = makedist('NegativeBinomial','R',r,'p',prob);
%set number of repetitions 
reps = 1000;
% set final time
T_final = 1.5;
T_treat=0.5;
T_stop=1;
% set up recording vectors
steps = T_final/15000;
t_vec=0:steps:T_final;[~,m]=size(t_vec);
A=zeros(m,reps);
time_vec=0:2*steps:T_final;

% loop over repetitions
for ii=1:reps
    mu_M=mu_M1;
    % set initial time
     t_s=0;
     % set recording step
    kk=1;
    % set initial populations size
    A_update=0;
    treat_num=0;
    while t_s<T_treat
        % if host is dead fill in recording vector with NaN and set time to
        % final
        if isnan(A_update)==1
            A(kk+1:end,ii)=NaN;
            t_s=T_final;
        else % if host is alive
          % set up propensity functions   
       propen1 = [psi; mu_M*A_update; mu_H+alpha*A_update];
       propen2=cumsum(propen1);
       prop0 =propen2(end);
       % generate a random number and use to find time jump
       r1 = rand(1);
       tau = 1/prop0*log(1/r1);
       % if time jump is smaller than a recording time step
       if tau<steps
           % update time
       t_s=t_s+tau;
       % generate random number to find process
       r2 = rand(1);
       propen1=1/prop0*propen1;
      
    if r2 < propen1(1)
        %if pick  up occurs generate clump size from dist
        C = random(pd1,1,1);
        % update current population
        A_update=A_update+C;
    elseif r2>propen1(1) && r2< propen1(2)+propen1(1)
        % if parasite death occurs reduce population by one
        A_update = A_update -1;
    else %if host dies set current value to NaN
        A_update = NaN;
    end
    % if new time is greater than next recording time update recording
    % vector
      if t_s> t_vec(kk+1)
          A(kk+1,ii)=A_update;
          kk=kk+1;
      end
       else % if time jump too big choose timestep length
           tau = steps;
           % update time
           t_s=t_s+tau;
           % no change in population
           A_update = A_update;
           % update recording vector
           A(kk+1,ii)=A_update;
          kk=kk+1;
       end
        end
   
    end
    if A_update>=C_num
        mu_M=mu_M1+u;
    end
    while t_s<T_stop
        
            % if host is dead fill in recording vector with NaN and set time to
        % final
        if isnan(A_update)==1
            A(kk+1:end,ii)=NaN;
            t_s=T_final;
        else % if host is alive
          % set up propensity functions   
       propen1 = [psi; mu_M*A_update; mu_H+alpha*A_update];
       propen2=cumsum(propen1);
       prop0 =propen2(end);
       % generate a random number and use to find time jump
       r1 = rand(1);
       tau = 1/prop0*log(1/r1);
       % if time jump is smaller than a recording time step
       if tau<steps
           % update time
       t_s=t_s+tau;
       % generate random number to find process
       r2 = rand(1);
       propen1=1/prop0*propen1;
      
    if r2 < propen1(1)
        %if pick  up occurs generate clump size from dist
        C = random(pd1,1,1);
        % update current population
        A_update=A_update+C;
    elseif r2>propen1(1) && r2< propen1(2)+propen1(1)
        % if parasite death occurs reduce population by one
        A_update = A_update -1;
    else %if host dies set current value to NaN
        A_update = NaN;
    end
    % if new time is greater than next recording time update recording
    % vector
      if t_s> t_vec(kk+1)
          A(kk+1,ii)=A_update;
          kk=kk+1;
      end
       else % if time jump too big choose timestep length
           tau = steps;
           % update time
           t_s=t_s+tau;
           % no change in population
           A_update = A_update;
           % update recording vector
           A(kk+1,ii)=A_update;
          kk=kk+1;
       end
        end
   
    end
    mu_M=mu_M1;
    
     while t_s<T_final
        
            % if host is dead fill in recording vector with NaN and set time to
        % final
        if isnan(A_update)==1
            A(kk+1:end,ii)=NaN;
            t_s=T_final;
        else % if host is alive
          % set up propensity functions   
       propen1 = [psi; mu_M*A_update; mu_H+alpha*A_update];
       propen2=cumsum(propen1);
       prop0 =propen2(end);
       % generate a random number and use to find time jump
       r1 = rand(1);
       tau = 1/prop0*log(1/r1);
       % if time jump is smaller than a recording time step
       if tau<steps
           % update time
       t_s=t_s+tau;
       % generate random number to find process
       r2 = rand(1);
       propen1=1/prop0*propen1;
      
    if r2 < propen1(1)
        %if pick  up occurs generate clump size from dist
        C = random(pd1,1,1);
        % update current population
        A_update=A_update+C;
    elseif r2>propen1(1) && r2< propen1(2)+propen1(1)
        % if parasite death occurs reduce population by one
        A_update = A_update -1;
    else %if host dies set current value to NaN
        A_update = NaN;
    end
    % if new time is greater than next recording time update recording
    % vector
      if t_s> t_vec(kk+1)
          A(kk+1,ii)=A_update;
          kk=kk+1;
      end
       else % if time jump too big choose timestep length
           tau = steps;
           % update time
           t_s=t_s+tau;
           % no change in population
           A_update = A_update;
           % update recording vector
           A(kk+1,ii)=A_update;
          kk=kk+1;
       end
        end
   
     end
end

 % find mean and variance of population size over time exluding NaNs for dead hosts
A_mean = nanmean(A,2);
A_var=nanvar(A,1,2);
t1=0:0.5/5000:0.5;
phi=psi;
% find pgf values for h(1), h'(1), h''(1) using h_neg_bin function
h1 = h_neg_bin(1);
%set up moment closure system 
 dx = @(t,x) [phi * (h1(2))  - alpha* x(2) - mu_M*x(1);phi*(h1(3)+ h1(2))+ mu_M*x(1)-2*mu_M*x(2)];
 % solve over time vector
 [~,x1]=ode45(dx,t1,[0;0]);
 
 M=x1(end,1);
 V=x1(end,2);
 % calculate p & r values
 p1=M/V;
 r1 = M^2/(V-M);
 % approximate distribution
 pdnew=makedist('NegativeBinomial','R',r1,'p',p1);
 % calculate truncated distributions at treatment level
 pdn2=truncate(pdnew,0,C_num);
 pdn3=truncate(pdnew,C_num,Inf);
 prop=cdf(pdnew,C_num);
 % find mean and variance of least and most infected from truncated
 % distributions
 M1=mean(pdn2);
 V1=var(pdn2);
 M2=mean(pdn3);
 V2=var(pdn3);
 t2=0.5:0.5/5000:1;
 % simulate untreated using mean and variance of least infected
 [~,x2]=ode45(dx,t2,[M1;V1]);
 % add in treatement
 mu_M=mu_M1+u;
 % simulate ode model for most infected
dx = @(t,x) [phi * (h1(2))  - alpha* x(2) - mu_M*x(1);phi*(h1(3)+ h1(2))+ mu_M*x(1)-2*mu_M*x(2)];
  [~,x3]=ode45(dx,t2,[M2;V2]);
  % calxculate combined mean and variance
  MT=prop*x2(:,1)+(1-prop)*x3(:,1);
  
  VT=(1-prop)*x2(:,2)+(prop)*x3(:,2)+prop*(1-prop)*(x2(:,1)-x2(:,1)).^2;

  x2=[MT,VT];
  mu_M=mu_M1;
  % simulate untretaed again
  dx = @(t,x) [phi * (h1(2))  - alpha* x(2) - mu_M*x(1);phi*(h1(3)+ h1(2))+ mu_M*x(1)-2*mu_M*x(2)];
 % solve over time vector
 t3=1:0.5/5000:1.5;
 [~,x3]=ode45(dx,t3,[MT(end);VT(end)]);
  t1=horzcat(t1,t2,t3);
 x1=vertcat(x1,x2,x3);
 % PLOT STOCHASTIC SIMULATION VERSUS ODE SIMULATIONS
 figure;
 subplot(2,1,1)
 plot(t_vec,A_mean,t1,x1(:,1),'LineWidth',4)
 xlabel('time')
 ylabel('mean')
 legend('simulation mean','approximated mean')
 subplot(2,1,2)
  plot(t_vec,A_var,t1,x1(:,2),'LineWidth',4)
   xlabel('time')
 ylabel('variance')
  legend('simulation variance','approximated variance')
  set(findall(gcf,'-property','FontSize'),'FontSize',15)
    
