%% code to simulate the Isham model using gillespie simulation and using
% deterministic moment closure model. Calculates error between two methods,
% creates plots and plots estimated distribution  based on negative
% binomial assumption.
close all
clear all
% set paramater values
global r p
psi = 52;
mu_H = 0;
mu_M= 10;
alpha=0.02;
% set mean clump size
M = 10;
%set tretament intensity
u=20;
%set clump variance
    r=0.5;    
    V=M.^2/r-M;
% set neg bin p value    
  prob=M/V;
  p=1-prob;
% set up clump distribution
pd1 = makedist('NegativeBinomial','R',r,'p',prob);
%set number of repetitions 
reps = 1000;
% set final time
T_final = 1;
% set up recording vectors
steps = T_final/5000;
t_vec=0:steps:T_final;
[~,m]=size(t_vec);
A=zeros(m,reps);
time_vec=0:2*steps:T_final;

% loop over repetitions
for ii=1:reps
    % set initial time
     t_s=0;
     % set recording step
    kk=1;
    % set initial populations size
    A_update=0;
    
  % when time is less than final time
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
err = zeros(2,5001);

% set up time vector
t = 0:1/5000:1;
%% Approximated model
phi=psi;
% find pgf values for h(1), h'(1), h''(1) using h_neg_bin function
h1 = h_neg_bin(1);
%set up moment closure system 
 dx = @(t,x) [phi * (h1(2))  - alpha* x(2) - mu_M*x(1); phi*(h1(3)+ h1(2))+ mu_M*x(1)-2*mu_M*x(2)];
 % solve over time vector
 [~,x1]=ode45(dx,t,[0;0]);
 err(1,:) = (A_mean-x1(:,1))./A_mean;
 err(2,:)=(A_var-x1(:,2))./A_var;
 %x1(:,2)=x1(:,2)-(x1(:,1).^2);
 %% plot code
% plot simulations v moment closure mean/variance
figure;
subplot(2,1,1)
plot(t_vec, A_mean, '-.r','LineWidth',4);
hold on
 plot(t(2:end),x1(2:end,1),'r','LineWidth',4);

xlabel('time')
ylabel('mean')
legend('simulation mean','moment closure mean')
subplot(2,1,2)
plot(t_vec, A_var, '-.b','LineWidth',4);
hold on
 plot(t(2:end),x1(2:end,2),'b','LineWidth',4);
xlabel('time')
ylabel('variance')
legend('simulation variance','moment closure variance')

figure;
subplot(2,1,1)
plot(t,err(1,:),'--r')
xlabel('time')
ylabel('relative error')
legend('relative error in mean')
subplot(2,1,2)
plot(t,err(1,:),'--b')
xlabel('time')
ylabel('relative error')
legend('relative error in variance')

% ste up probability dist with p and r values determined by moment
% closure system mean and variance
p2=x1(end,1)/x1(end,2);
r2=x1(end,1)^2/(x1(end,2)-x1(end,1));



pd2=makedist('NegativeBinomial','R',r2,'p',p2);

% take large sample from prob dist
sam2=random(pd2,1,10000);
% take histogram data of distribution of simulations at two separate time
% points later enough to be at equilibrium
[c1,v1]=hist(A(3500,:),0:5:160);
[c3,v3]=hist(A(1500,:),0:5:160);
% take hisogram of moment closure sample
[c2,v2]=hist(sam2,0:5:160);
% turn histogram counts into proportions
c1=c1/sum(c1);
c3=c3/sum(c3);
c2=c2/sum(c2);
v1=[0,v1];
v2=[0,v2];
c1=[0,c1];
c2=[0,c2];
v3=[0,v3];
c3=[0,c3];
% plot moment closure distribution v actual distribution from simulation
% sample times

figure;
plot(v1,c1,'r',v3,c3,'g',v2,c2,'b','LineWidth',4)
xlabel('parasite burden')
ylabel('proportion of hosts burdened')
legend('exact stochastic 1 ','exact stochastic 2', 'moment closure approximation')

%% seection of 10 percent of hosts
B=A(1500,:);
B=B(~isnan(B));
m=length(B);
pr=round(m/10);
%select at random
r1=randi([1, m],[1,pr]);

C=B(r1);
%select most infected
B=sort(B,'descend');
D=B(1:pr);

%plot distribution in these hosts
[c1,v1]=hist(C,0:10:200);
[c2,v2]=hist(D,0:10:200);
c1=c1/sum(c1);
c2=c2/sum(c2);
figure;
plot(v1,c1,'r',v2,c2,'b','LineWidth',4)
xlabel('number of parasites')
ylabel('proportion of hosts')
legend('10%randomly selected','Most infected 10% selected')

