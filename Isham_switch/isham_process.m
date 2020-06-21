clear all
close all
global phi alpha mu_M p r u tu
% set parameter values
phi = 52/12;
mu_H = 0;
mu_M= 1/12;
u=20;
alpha=0.002;
% set mean clump size
M = 0.5;
% set dispersion parameter
dx=@(t,x)isham(t,x);
    r=0.5;    
% fin variance in clump sizes
init=[25,45];
   
V = (M.^2+r*M)/r;
   %find p parameter 
  prob=M/V;
  p=1-prob;
  % set treatment time vector
tu = 0:0.25:3.5;
m=length(tu);
% find all combos of treat/no treat in time intervals
u1 = de2bi(0:2^(m)-1);
% set cost function parameters
alp1=0/(2*11.51^2);
alp2=2/(2*20.755^2);
MT=1.98;
VT=3.49;
[n,~]=size(u1);
% set initial cost to infinity for comparison
costt1=inf;
% loop over all treatment combos
for ii=1:n
    
    u=u1(ii,:);
    ts=0:0.005:3.5;
    % simulate chosen combo on model
    [~,xs]=ode45(dx,ts,init);
    % calculate cost
    Y1=alp1*(xs(:,1)-MT).^2+alp2*(xs(:,2)-VT).^2;
    cost=0.005*trapz(Y1)+0.25*sum(u(1:end-1));
    % check if cost less than current best if yes replace current best
    if cost < costt1
        t_num1=ii;
        costt1=cost;
         t1=ts;
        x1=xs;
    end
       
    
    
end
% plot final best option and response to control
costu1=0.25*sum(u1(t_num1,1:end-2));
figure;
subplot(3,1,1)
stairs(tu,u1(t_num1,:),'k','LineWidth',4)
xlabel('time')
ylabel('control on/off')
subplot(3,1,2)
plot(t1,x1(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('mean')
subplot(3,1,3)
plot(t1,x1(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
%% repeat with alternative time vector

tu = 0:0.5:3.5;
m=length(tu);
u2 = de2bi(0:2^(m)-1);
[n,~]=size(u2);
costt2=inf;
for ii=1:n
    
    u=u2(ii,:);
    ts=0:0.005:3.5;
    [~,xs]=ode45(dx,ts,init);
    Y1=alp1*(xs(:,1)-MT).^2+alp2*(xs(:,2)-VT).^2;
    cost=0.005*trapz(Y1)+0.5*sum(u(1:end-1));
    if cost < costt2
        t_num2=ii;
        costt2=cost;
         t1=ts;
        x1=xs;
    end
       
    
    
end
costu2=0.5*sum(u2(t_num2,1:end-2));
figure;
subplot(3,1,1)
stairs(tu,u2(t_num2,:),'k','LineWidth',4)
xlabel('time')
ylabel('control on/off')
subplot(3,1,2)
plot(t1,x1(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('mean')
subplot(3,1,3)
plot(t1,x1(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
%% repeat with alternative timevector
tu = 0:1:4;
m=length(tu);
u3 = de2bi(0:2^(m)-1);
[n,~]=size(u3);
costt3=inf;
for ii=1:n
    
    u=u3(ii,:);
    ts=0:0.005:3.5;
    [~,xs]=ode45(dx,ts,init);
    Y1=alp1*(xs(:,1)-MT).^2+alp2*(xs(:,2)-VT).^2;
    cost=0.005*trapz(Y1)+sum(u(1:end-1))-0.5*u(end-1);
    if cost < costt3
        t_num3=ii;
        costt3=cost;
         t1=ts;
        x1=xs;
    end
       
    
    
end
costu3=sum(u3(t_num3,1:end-2))+0.5*u3(t_num3,end-1);
figure;
subplot(3,1,1)
stairs(tu,u3(t_num3,:),'k','LineWidth',4)
xlabel('time')
ylabel('control on/off')
subplot(3,1,2)
plot(t1,x1(:,1),'r','LineWidth',4)
xlabel('time')
ylabel('mean')
subplot(3,1,3)
plot(t1,x1(:,2),'b','LineWidth',4)
xlabel('time')
ylabel('variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

%% sequential comparison on same time vectors

t1=0:0.25:3.5;
u1=zeros(size(t1));
m=length(t1);
init1= init;
costus1=0;
costr1=0;
for ii=1:m-1
    u=[0,0];
    tu=[t1(ii), t1(ii+1)];
    
        [ts1,xs1]=ode45(dx,[t1(ii),t1(ii+1)],init1);
        u=[1 1];
         [ts2,xs2]=ode45(dx,[t1(ii),t1(ii+1)],init1);
         
         cost1=alp1*(xs1(:,1)-MT).^2+alp2*(xs1(:,2)-VT).^2;
         cost1=trapz(cost1,ts1);
         
         cost2=alp1*(xs2(:,1)-MT).^2+alp2*(xs2(:,2)-VT).^2;
         
        cost2=trapz(cost2,ts2)+0.25;
       if cost1<cost2
           u1(ii)=0;
           init1=xs1(end,:);
          
           
       else
           u1(ii)=1;
           init1=xs2(end,:);
          
           
       end
end
costus1=0.25*sum(u1(1:end-1));
  tu=t1;
u=u1;
t1=0:0.05:3.5;
[ts,xs]=ode45(dx,t1,init);
cost2=alp1*(xs(:,1)-MT).^2+alp2*(xs(:,2)-VT).^2;
costr1=0.05*trapz(cost2);

t1=0:0.5:3.5;
u1=zeros(size(t1));
m=length(t1);
init1= init;
costus2=0;
costr2=0;
for ii=1:m-1
    u=[0,0];
    tu=[t1(ii), t1(ii+1)];
        [ts1,xs1]=ode45(dx,[t1(ii),t1(ii+1)],init1);
        u=[1 1];
         [ts2,xs2]=ode45(dx,[t1(ii),t1(ii+1)],init1);
         cost1=alp1*(xs1(:,1)-MT).^2+alp2*(xs1(:,2)-VT).^2;
         cost1=trapz(cost1,ts1);
         
         cost2=alp1*(xs2(:,1)-MT).^2+alp2*(xs2(:,2)-VT).^2;
         
        cost2=trapz(cost2,ts2)+(ts2(end)-ts2(1));
       if cost1<cost2
           u1(ii)=0;
           init1=xs1(end,:);
      
       else
           u1(ii)=1;
           init1=xs2(end,:);
                
       end
end
costus2=0.5*sum(u1(1:end-1));
tu=t1;
u=u1;
t1=0:0.05:3.5;
[ts,xs]=ode45(dx,t1,init);
cost2=alp1*(xs(:,1)-MT).^2+alp2*(xs(:,2)-VT).^2;
costr2=0.05*trapz(cost2);
t1=0:1:4;
u1=zeros(size(t1));
m=length(t1);
init1= init;
costus3=0;
costr3=0;
for ii=1:m-1
    if ii<m-1
    u=[0,0];
    tu=[t1(ii), t1(ii+1)];
        [ts1,xs1]=ode45(dx,[t1(ii),t1(ii+1)],init1);
        u=[1 1];
         [ts2,xs2]=ode45(dx,[t1(ii),t1(ii+1)],init1);
    else
         u=[0,0];
    tu=[t1(ii), t1(ii)+0.5];
        [ts1,xs1]=ode45(dx,[t1(ii),t1(ii)+0.5],init1);
        u=[1 1];
         [ts2,xs2]=ode45(dx,[t1(ii),t1(ii)+0.5],init1);
    end
        cost1=alp1*(xs1(:,1)-MT).^2+alp2*(xs1(:,2)-VT).^2;
         cost1=trapz(cost1,ts1);
         
         cost2=alp1*(xs2(:,1)-MT).^2+alp2*(xs2(:,2)-VT).^2;
        cost2=trapz(cost2,ts2)+1*(ts2(end)-ts2(1));
       if cost1<cost2
           u1(ii)=0;
           init1=xs1(end,:);
         
       else
           u1(ii)=1;
           init1=xs2(end,:);
          
          
       end
end
costus3=sum(u1(1:end-1))-0.5*(u1(end-1));
tu=t1;
u=u1;
t1=0:0.05:3.5;
[ts,xs]=ode45(dx,t1,init);
cost2=alp1*(xs(:,1)-MT).^2+alp2*(xs(:,2)-VT).^2;
costr3=0.05*trapz(cost2);


%% create cost vector and plot as bar chart for comparisons between treatment interval lengths and optimal vs sequential
figure;
cost_v=[costt3-costu3,costu3;costr3,costus3;nan,nan;costt2-costu2,costu2;costr2,costus2;nan,nan;costt1-costu1,costu1;costr1,costus1];
bar(cost_v,'stacked');
hold on;
ax = gca;
ax.XTickLabels = {'\tau=1','\tau=1, seq','','\tau=0.5','\tau=0.5, seq','','\tau=0.25','\tau=0.25, seq'};
grid on;
xlabel('Discrete Interval length')
ylabel('Cost')
set(findall(gcf,'-property','FontSize'),'FontSize',15)