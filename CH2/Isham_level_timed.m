%% performs simulations of Isham model with treatment which begins when a 
% host has burden exceeding value C_max and ends after time given by tau1.
% Creates solution plots and compares with results of untreated and fully
% treated models.
% set parameters
global r p
psi = 52;
mu_H = 0;
mu_M= 10;
alpha=0.02;
M = 10;
% set death due to treatment rate
u=20;
% set time length of a single treatment
tau1 = 0.05;
% set low enough value that treatment is off at start for time last
% treatment started
t_star=-0.5;
% set trigger value for treatment
Cmax = 50;
% set up probability distribution fo rclump  sizes
M = 10;

    r=0.5;    


    V=M.^2/r-M;
    
  prob=M/V;
  p=1-prob;

pd1 = makedist('NegativeBinomial','R',r,'p',prob);
%set number of repetitions
reps = 500;
%ste final time
T_final = 1;
% set up recording vectors
steps = T_final/5000;
t_vec=0:steps:T_final;
[~,m]=size(t_vec);
A=zeros(m,reps);
%loop over repetitions
for ii=1:reps
    %set initial time and populations
    
    t_s=0;
    kk=1;
    A_update=0;
    %set up while loop for time 
    while t_s<T_final
        % if host is dead pop= NaN fill recording vector and ste time to
        % final time
        if isnan(A_update)==1
            A(kk+1:end,ii)=NaN;
            t_s=T_final;
        else %if host is not dead
       if A_update >=Cmax %if population exceeds trigger value
           %record treatment start time
           t_star=t_s;
           % set up propensity functions increasd death rate
            propen1 = [psi; (mu_M+u)*A_update; mu_H+alpha*A_update];
           
       elseif t_s<=t_star+tau1 %if population is already being treated
           % set up propensity functions with increasd death rate
           propen1 = [psi; (mu_M+u)*A_update; mu_H+alpha*A_update];
           
       else %if host not being treated
           
         %set up propensity functions with standard parasite death rate  
       propen1 = [psi; mu_M*A_update; mu_H+alpha*A_update];
       end
       propen2=cumsum(propen1);
       prop0 =propen2(end);
       %generate random number an dgenerate timeuntil next event
       r1 = rand(1);
       tau = 1/prop0*log(1/r1);
       if tau < steps % if time step shorter than recording vector step length
           %update time
       t_s=t_s+tau;
       %find what event happens
       r2 = rand(1);
       propen1=1/prop0*propen1;
    if r2 < propen1(1)
        %if pick up generate clump size anfd update current populations
        %size
        C = random(pd1,1,1);
        A_update=A_update+C;
    elseif r2>propen1(1) && r2< propen1(2)+propen1(1)
        %if parasite death update population by removing one
        A_update = A_update -1;
    else
        %if host death occurs set population size to NaN
        A_update = NaN;
    end
      if t_s> t_vec(kk+1)
          %if current time greater than next recording vector time updatae
          %recording vector with population size
          A(kk+1,ii)=A_update;
          kk=kk+1;
      end
       else
           %if time step too big set to lower value 
           tau = steps;
           %update time
           t_s=t_s+tau;
           %update recording vectors
           A_update = A_update;
           A(kk+1,ii)=A_update;
          kk=kk+1;
       end
        end
    end
    
end
% find mean and variance of simulations excluding dead hosts
A_mean = nanmean(A,2);
A_var=nanvar(A,1,2);

t = 0:0.01:1;
% find moment closure solution for mean and variuance
phi=psi;
h1 = h_neg_bin(1);
 dx = @(t,x) [phi * (h1(2))  - alpha* x(2) - mu_M*x(1); phi*(h1(3)+ h1(2))+ mu_M*x(1)-2*mu_M*x(2)];
 
 [~,x1]=ode45(dx,t,[0;0]);
 dx2 = @(t,x) [phi * (h1(2))  - alpha* x(2) - (mu_M+u)*x(1); phi*(h1(3)+ h1(2))+ (mu_M+u)*x(1)-2*(mu_M+u)*x(2)];
 
 [~,x2]=ode45(dx2,t,[0;0]);
 
%% distributions
% set up distriutions for moment closure estimates
%untreated dist
p2=x1(end,1)/x1(end,2);
r2=x1(end,1)^2/(x1(end,2)-x1(end,1));

%treated dist
p3=x2(end,1)/x2(end,2);
r3=x2(end,1)^2/(x2(end,2)-x2(end,1));


pd3=makedist('NegativeBinomial','R',r3,'p',p3);

pd2=makedist('NegativeBinomial','R',r2,'p',p2);
% sample from estimate distributions
sam1=random(pd3,1, 10000);
sam2=random(pd2,1,10000);
% histogram samples and examples of simulations at time couple of time
% points
% sampled distributions from simulations
[c1,v1]=hist(A(end,:),0:5:160);
[c3,v3]=hist(A(2000,:),0:5:160);

%distrbutions sampled from approximated distributions
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
