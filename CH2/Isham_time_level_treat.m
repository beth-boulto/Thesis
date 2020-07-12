%% code to simulate Isham model with treatment which occurs during set
% intervals if the host burden exceeds C_max. Simulated usign a gillespie
% simulation with plots to ompare with treated and untread moment closure
% approximations.
close all
clear all
% set up parameters
global r p
psi = 25;
mu_H = 0;
mu_M1=10;
u=20;
alpha=0.02;
M = 10;
%set burden limit
Cmax=50;
%set clump dist parameters
    r=0.5;    
    V=M.^2/r-M;
     prob=M/V;
  p=1-prob;

pd1 = makedist('NegativeBinomial','R',r,'p',prob);
%set repetitions
reps = 500;
%set final time
T_final = 1;
% set up recording vectors
steps = T_final/5000;
t_vec=0:steps:T_final;
[~,m]=size(t_vec);
A=zeros(m,reps);
% set up when differnt treatment periods occur
Time_change = 0:2*steps:T_final;
[~,m]=size(Time_change);
% loop over repetitions
for ii=1:reps
    % initialise time and population
     t_s=0;
    kk=1;
    A_update=0;
    %loop over treatment periods
    for jj=2:m
        %set parasite death rate for current treatment period
        if A_update>Cmax
            % treated
            mu_M=mu_M1+u;
        else
            %untreated
            mu_M=mu_M1;
        end
        while t_s<Time_change(jj) % run until next treatment period
            if isnan(A_update)==1 % if host is dead end simulation record NaNs for population size
            A(kk+1:end,ii)=NaN;
            t_s=T_final;
            else
                % set up propensity functions
             propen1 = [psi; mu_M*A_update; mu_H+alpha*A_update];
             propen2=cumsum(propen1);
       prop0 =propen2(end);
       % random number generation to find time untilnext event
       r1 = rand(1);
       tau = 1/prop0*log(1/r1);
       tp = t_vec(kk+1)-t_s;
       if tau<tp % if time jump smaller than time to next recording step
           %update time
           t_s=t_s+tau;
           % generarte random number to find what process happens
       r2 = rand(1);
       propen1=1/prop0*propen1;
    if r2 < propen1(1)
        %if birth find clump size and update
        C = random(pd1,1,1);
        A_update=A_update+C;
    elseif r2>propen1(1) && r2< propen1(2)+propen1(1)
        %if parasite death remove parasite from pop
        A_update = A_update -1;
    else
   % if host dies set popo to NaN
        A_update = NaN;
    end
       else
           %if time step too great choose smaller step and update time and
           %population
           tau=tp;
           t_s=t_s+tau;
           A(kk+1,ii)=A_update;
          kk=kk+1;
       end
           
       
            end
        end
    end
    
end
% find simulation mean and variance excluding NaNs
 A_mean = nanmean(A,2);
A_var=nanvar(A,1,2);

t = 0:0.01:1;
% moment closure solution of untreated and fully treated simulation modoels
phi=psi;
h1 = h_neg_bin(1);
 dx = @(t,x) [phi * (h1(2))  - alpha* x(2) - mu_M*x(1); phi*(h1(3)+ h1(2))+ mu_M*x(1)-2*mu_M*x(2)];
 
 [~,x1]=ode45(dx,t,[0;0]);
 dx2 = @(t,x) [phi * (h1(2))  - alpha* x(2) - (mu_M+u)*x(1); phi*(h1(3)+ h1(2))+ (mu_M+u)*x(1)-2*(mu_M+u)*x(2)];
 
 [~,x2]=ode45(dx2,t,[0;0]);
 

% set up distriutions fo rmoment closure estimates
p2=x1(end,1)/x1(end,2);
r2=x1(end,1)^2/(x1(end,2)-x1(end,1));

p3=x2(end,1)/x2(end,2);
r3=x2(end,1)^2/(x2(end,2)-x2(end,1));


pd3=makedist('NegativeBinomial','R',r3,'p',p3);

pd2=makedist('NegativeBinomial','R',r2,'p',p2);
% sample from estimate distributions
sam1=random(pd3,1, 10000);
sam2=random(pd2,1,10000);
% histogram samples and examples of simulations at time couple of time
% points
[c1,v1]=hist(A(end,:),0:5:160);
[c3,v3]=hist(A(2000,:),0:5:160);
[c2,v2]=hist(sam2,0:5:160);
[c4,v4]=hist(sam1,0:5:160);
% set histogram counts to proportions
c1=c1/sum(c1);
c3=c3/sum(c3);
c2=c2/sum(c2);
c4=c4/sum(c4);

v1=[0,v1];
v2=[0,v2];
c1=[0,c1];
c2=[0,c2];
v3=[0,v3];
c3=[0,c3];
v4=[0,v4];
c4=[0,c4];

