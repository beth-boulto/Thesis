%% code which performs gillespie simulation to model the larval model from
%the end of chapter 2 and plot results and comparisons with expected
%distribution based on negative binomial distribution with parameters
%calculated form mean and variance.
clear all
close all

%set parameters
N = 100;
mu_L=2;
mu_M=1;
phi=03/N;
beta =300;                                                                                                                                                                                                                                                             
reps = 50;
M = 2;
%set mean clump size for pick up
D=2;
D_M=1/20;
r2=0.1; 
    r=0.5;    
%set clump pick up variance
V = (M.^2+r*M)/r;
  prob=M/V;
  p=1-prob;
% set up clumped pick up distribution
pd1 = makedist('NegativeBinomial','R',r,'p',prob);
% set up recording vectors
T_final=10;
steps=5000;
t_vec=0:1/steps:T_final;
[~,m]=size(t_vec);

A=zeros(N+1,m);
%set parasite burdens to initial value of zero
Init = zeros(N+1,1);
% set larval pool to 500 initially
Init(1)=500;
% repeat over number of simulatiojns
for ii=1:reps
    % set up initial conditions time, stated etc
    A(:,1,ii)=Init;
    t_s=0;
    A_update=Init;
    kk=1;
    N_t=N;
    index=1:N;
    % loop over time less than end point
    while t_s<T_final
        % set up propensity functions
       propen1 = [N*phi*A_update(1);mu_L*A_update(1);beta/N*sum(A_update(2:end))];
       propen2 = mu_M*A_update(2:end)+D_M*A_update(2:end).^2;
       propen1=[propen1;propen2];
       prop1 = cumsum(propen1);
       prop0=prop1(end);
       % generate random number and find time jump
       r1 = rand(1);
       tau = 1/prop0*log(1/r1);
       if prop0 ==0
           % if no parasites remain jump to end
           A(:,kk+1:m) = A_update;
           t_s=T_final;
       else
          if tau<1/steps 
      % if time jump small enough move time forward
       t_s=t_s+tau;
         
       % find which event occurs   
       r2 = rand(1);
       prop1=1/prop0*prop1;
       proc=find(r2 < prop1, 1, 'first');
       if proc==1
           % if pick up find number of parasites picked up and which host
           % gets them
           pd2 = truncate(pd1,0,A_update(1));
           C=random(pd2);
           r3 = randi(N_t);
           A_update(r3+1)=A_update(r3+1)+C;
           A_update(1)=A_update(1)-C;
       elseif proc==2
           % if larval death kill off one larva
           A_update(1)=A_update(1)-1;
       elseif proc==3 
           % if birth event add more larvae
           A_update(1)=A_update(1)+D;
             
       else
           %if mature parasite death kill off parasite in host 
           A_update(proc-2)=A_update(proc-2)-1;
       end
        elseif tau<0
              return
          else
           % if time jump too large choose smaller leap update time
           tau = 1/steps;
           t_s = t_s+tau;
         
          end
            if t_s>t_vec(kk+1)
                % if time to record update recording vectors
           A(:,kk+1,ii)=A_update;
           kk=kk+1;
            end
       
       end
    end
   
   % find mean and variance of mature parasite populations 
   B=A(2:end,:,ii);
    A_mean(ii,:) = nanmean(B,1);
    A_var(ii,:)=nanvar(B,1);
    % end loop over simulations
end
% find mean and variance in larval population over simulations
larv1 = zeros(reps,m);
for ii=1:reps
    larv1(ii,:)=A(1,:,ii);
    larv1(ii,1)=500;
end
larv=zeros(2,m);
%calculate mean and variance of each simulation
larv(1,:)=mean(larv1,1);
larv(2,:)=var(larv1,1);

% plot larval mean and variance over simulations
figure;
plot(t_vec,larv(1,:),'r',t_vec,larv(2,:),'b','LineWidth',4)
xlabel('time')
ylabel('larval mean/variance')
legend('larval mean', 'larval variance')
% find mean of mature mean and variance of mature mean over simulations
A_mm = mean(A_mean,1);
A_mv=var(A_mean,1);

% plot mean and variance over simulations 
figure;
plot(t_vec,A_mm,'r',t_vec,A_mv,'b','LineWidth',4)
xlabel('time')
ylabel('mean/variance')
legend('mean of mature means','variance of mature means')

plot(t_vec,larv(1,:),'g','LineWidth',4)
xlabel('time')
ylabel('mean')
legend('larval mean')
%% code for plotting indiviual simulations
figure;
plot(t_vec,A(1,:,50),'g','LineWidth',4,'LineWidth',4)
xlabel('time')
ylabel('larval simulation')

figure;
plot(t_vec,A_mean(50,:),'r',t_vec,A_var(50,:),'b','LineWidth',4)
xlabel('time')
ylabel('mean and variance')
legend('mean burden', 'variance in burden')
% distributions of individual simulations as histogram plots
[counts,val]=hist(A(2:end,25000,10),0:5:150);
[c2,v2]=hist(A(2:end,30000,10),0:5:150);
[c3,v3]=hist(A(2:end,end,10),0:5:150);
%normalise to give proportions
counts=counts/(sum(counts));
c2=c2/sum(c2);
c3=c3/sum(c3);
%plot distributions
figure;
plot(0:2:150,counts,'r',0:2:150,c2,'b',0:2:150,c3,'g','LineWidth',4)
xlabel('parasite burden')
ylabel('proportion of hosts')



%% new plot codes - mean and variances vs mean of mean and var

figure;

hold on
for ii=1:reps
    plot(t_vec,larv1(ii,:))
end
    plot(t_vec,larv(1,:),'g','LineWidth',4) 
xlabel('time')
ylabel('larval population')


figure
hold on
for ii=1:reps
    plot(t_vec,A_mean(ii,:))
end
plot(t_vec,A_mm,'r','LineWidth',4) 
xlabel('time')
ylabel('distribution mean')





figure;
hold on
for ii=1:reps
    plot(t_vec,A_var(ii,:))
end
plot(t_vec,A_mv,'--b','LineWidth',4) 
xlabel('time')
ylabel('distribution variance')


p=zeros(reps,m);
r=zeros(reps,m);
%calculate p and r of individual simulation at each time recorded point
for ii=1:reps
    p(ii,:)=A_mean(ii,:)./A_var(ii,:);
    r(ii,:)=A_mean(ii,:).^2./(A_var(ii,:)-A_mean(ii,:));
end

%calculate mean p and r values ( pm2 & rm2 are correct versions
pm=A_mm./A_vm;
rm=A_mm.^2./(A_vm-A_mm);

pm2 =A_mm./A_mv;
rm2 = A_mm.^2./(A_mv-A_mm);
%plot figures for p and r values of simulations vs means  
figure;

hold on
for ii=1:reps
    plot(t_vec,p(ii,:))
end
plot(t_vec,pm,'--k','LineWidth',4) 
xlabel('time')
ylabel('p parameter')

figure;

hold on
for ii=1:reps
    plot(t_vec,r(ii,:))
end
plot(t_vec,rm,'--k','LineWidth',4) 
xlabel('time')
ylabel('r parameter')

figure;

hold on
for ii=1:reps
    plot(t_vec,p(ii,:))
end
plot(t_vec,pm2,'--k','LineWidth',4) 
xlabel('time')
ylabel('p parameter')

figure;

hold on
for ii=1:reps
    plot(t_vec,r(ii,:))
end
plot(t_vec,rm2,'--k','LineWidth',4) 
xlabel('time')
ylabel('r parameter')





%% average dist
%find p & r values of neg bin dist
p=pm(end);
r=rm(end);

pd1=makedist('NegativeBinomial','R',r,'p',p);

% sample from distribution
sam1 = random(pd1,[1 10000]);

%choose mean and variaince of 10th run of simulation at the end of the time
m=A_mean(10,end);
v=A_var(10,end);
%calculate p & r for 10th simualtion at end time
p=m/v
r=m^2/(v-m)

pd1=makedist('NegativeBinomial','R',r,'p',p);
%sample from resulting distribution
sam2 = random(pd1,[1 10000]);
%repeat or 30th simulation
m=A_mean(30,end);
v=A_var(30,end);
p=m/v;
r=m^2/(v-m);

pd1=makedist('NegativeBinomial','R',r,'p',p);

sam3 = random(pd1,[1 10000]);
%calculate distribution o samples
[c1,v1]=hist(sam1,0:5:160);
[c2,v2]=hist(sam2,0:5:160);
[c3,v3]=hist(sam3,0:5:160);

c1=c1/sum(c1);
c3=c3/sum(c3);
c2=c2/sum(c2);
v1=[0,v1];
v2=[0,v2];
c1=[0,c1];
c2=[0,c2];
v3=[0,v3];
c3=[0,c3];

%plot sample distributions
figure;
plot(v1,c1,'k',v2,c2,'r',v3,c3,'b','LineWidth',4)
xlabel('number of parasites')
ylabel('proportion of hosts')
legend('average distribution', 'example distribution 1','example distribution 2')


