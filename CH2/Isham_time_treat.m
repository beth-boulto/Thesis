%% code to simulate the Isham model stochasticallly and using the moment
%clousre approximation with a treatment of all hosts that begin halfway
%through the overall time period and ends at the end of the time period.
%Plots comparisons of two methods and expected distributions based on
%negative binomialdistribution.
close all
clear all
global r p
% set paramaeters
psi = 52;
mu_H = 0;
mu_M= 10;
alpha=0.02;
% set clump mean and variance
M = 10;
u=20;
    r=0.5;    


    V=M.^2/r-M;
    
  prob=M/V;
  p=1-prob;
% set clump distribution
pd1 = makedist('NegativeBinomial','R',r,'p',prob);
% set repetrions
reps = 100;
% set time to start treatment and end time
T_final1 = 0.5;
T_final2=1;
% set time vector
steps = T_final2/5000;
t_vec=0:steps:T_final2;
[~,m]=size(t_vec);
% ste recording vector
A=zeros(m,reps);
%start simulations
for ii=1:reps
    t_s=0;
    kk=1;
    A_update=0;
    while t_s<T_final1
        % check for living host
        if isnan(A_update)==1
            % if dead fill remaining recording vector with nan and end rep
            A(kk+1:end,ii)=NaN;
            t_s=T_final1;
        else
            %otherwise- set up propensity functios
       propen1 = [psi; mu_M*A_update; mu_H+alpha*A_update];
       propen2=cumsum(propen1);
       % set p_0
       prop0 =propen2(end);
       % simulate time jump
       r1 = rand(1);
       tau = 1/prop0*log(1/r1);
       if tau < steps
       t_s=t_s+tau;
       % determine event which occurs
       r2 = rand(1);
       propen1=1/prop0*propen1;
    if r2 < propen1(1)
        %pick up
        C = random(pd1,1,1);
        A_update=A_update+C;
    elseif r2>propen1(1) && r2< propen1(2)+propen1(1)
        %parasite death
        A_update = A_update -1;
    else
        %host death
        A_update = NaN;
    end
      if t_s> t_vec(kk+1)
          %update recording vector
          A(kk+1,ii)=A_update;
          kk=kk+1;
      end
       else
           tau = steps;
           t_s=t_s+tau;
           A_update = A_update;
           A(kk+1,ii)=A_update;
          kk=kk+1;
       end
        end
    end
    %end of untreated section
    %begin treated sections - repeat of above but mu_M = mu_M+u now for all
    %hosts
     while t_s<T_final2 && t_s>=T_final1
        if isnan(A_update)==1
            A(kk+1:end,ii)=NaN;
            t_s=T_final2;
        else
       %calculate propensity functions     
       propen1 = [psi; (mu_M+u)*A_update; mu_H+alpha*A_update];
       propen2=cumsum(propen1);
       prop0 =propen2(end);
       %generate random number and calcute time jump
       r1 = rand(1);
       tau = 1/prop0*log(1/r1);
       %check for time step length and recording times
       if tau < steps
       t_s=t_s+tau;
       %simulate which event has occurred
       r2 = rand(1);
       propen1=1/prop0*propen1;
    if r2 < propen1(1)
        %clumped pick up update
        C = random(pd1,1,1);
        A_update=A_update+C;
    elseif r2>propen1(1) && r2< propen1(2)+propen1(1)
        %death update
        A_update = A_update -1;
    else
        %host dies
        A_update = NaN;
    end
      if t_s> t_vec(kk+1)
          A(kk+1,ii)=A_update;
          kk=kk+1;
      end
       else
           tau = steps;
           t_s=t_s+tau;
           A_update = A_update;
           A(kk+1,ii)=A_update;
          kk=kk+1;
       end
        end
    end
    
end
%calculate mean and variance excluding dead hosts
A_mean = nanmean(A,2);
A_var=nanvar(A,1,2);

% moment closure approximated ode
phi=psi;
t = 0:0.01:0.5;
t2 = 0.5:0.01:1;
% set pgf of neg bin and derivatives
h1 = h_neg_bin(1);
% set untreated ode model
 dx = @(t,x) [phi * (h1(2))  - alpha* x(2) - mu_M*x(1); phi*(h1(3)+ h1(2))+ mu_M*x(1)-2*mu_M*x(2)];
 % simualte untreated
 [~,x1]=ode45(dx,t,[0;0]);
 % set up treated ode model
 dx = @(t,x) [phi * (h1(2))  - alpha* x(2) - (mu_M+u)*x(1); phi*(h1(3)+ h1(2))+ (mu_M+u)*x(1)-2*(mu_M+u)*x(2)];
 % calculate approximate p and r values of untreated distribution
 p2=x1(end,1)/x1(end,2);
r2=x1(end,1)^2/(x1(end,2)-x1(end,1));
% simulate tretaed 
 [~,x2]=ode45(dx,t2,x1(end,:));
 t=[t,t2];
 x1=vertcat(x1,x2);

% calculate p & r of stochastic simulations
A_m = mean(A_mean(3500:end));
A_v=mean(A_var(3500:end));
p2=x1(50,1)/x1(50,2);
r2=x1(50,1)^2/(x1(50,2)-x1(50,1));


% set up approximate distributions from p/r
pd2=makedist('NegativeBinomial','R',r2,'p',p2);

% smapke from distribution and histogram
sam1=random(pd1,1, 10000);
sam2=random(pd2,1,10000);
[c1,v1]=hist(A(end,:),0:5:160);
[c3,v3]=hist(A(4000,:),0:5:160);
[c2,v2]=hist(sam2,0:5:160);
c1=c1/sum(c1);
c3=c3/sum(c3);
c2=c2/sum(c2);
sum(c1)
sum(c3)
v1=[0,v1];
v2=[0,v2];
c1=[0,c1];
c2=[0,c2];
v3=[0,v3];
c3=[0,c3];
% set up ode distribtiotn from p/r
p2=x1(50,1)/x1(50,2);
r2=x1(50,1)^2/(x1(50,2)-x1(50,1));
pd2=makedist('NegativeBinomial','R',r2,'p',p2);
%sample and histogram
sam1=random(pd1,1, 10000);
[c4,v4]=hist(sam1,0:5:160);
c4=c4/sum(c4);
v4=[0,v4];
c4=[0,c4];
pd2=makedist('NegativeBinomial','R',r2,'p',p2);
