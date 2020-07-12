% code for examining phase plane portrait of variance inclusive model with
% and without the application of treatment
clear all 
close all
global beta mu_M D_M MT VT T_start
beta=1.5;
D_M=1/20;
mu_M1=1/(5*365);

mu_M=mu_M1;

dx=@(t,x)metapop(t,x);
%% phase plane no treatment
[X,Y]=meshgrid(0.5:1:30.5,0.1:0.5:15.1);
u=zeros(size(X));
v=zeros(size(X));
u1=zeros(size(X));
v1=zeros(size(X));
for ii=1:numel(X)
    dv=dx(1,[X(ii),Y(ii)]);
    u(ii)=dv(1)./((dv(1)^2+dv(2)^2)^(1/2));
    v(ii)=dv(2)./((dv(1)^2+dv(2)^2)^(1/2));
    u1(ii)=dv(1)/X(ii);
    v1(ii)=dv(2)/Y(ii);
end

[ts,xs]=ode45(dx,[0 5],[9;7]);
figure;
hold on
quiver(X,Y,u,v,'r')
plot(xs(:,1),xs(:,2),'g','Linewidth',4)
xlabel('mean')
ylabel('variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
figure;
quiver(X,Y,u1,v1,'b')

xlabel('mean')
ylabel('variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
%% phase plane during treatment
mu_M=mu_M1+1;
[X,Y]=meshgrid(0.5:1:30.5,0.1:0.5:15.1);
u=zeros(size(X));
v=zeros(size(X));
u1=zeros(size(X));
v1=zeros(size(X));
for ii=1:numel(X)
    dv=dx(1,[X(ii),Y(ii)]);
    u(ii)=dv(1)./((dv(1)^2+dv(2)^2)^(1/2));
    v(ii)=dv(2)./((dv(1)^2+dv(2)^2)^(1/2));
    u1(ii)=dv(1)/X(ii);
    v1(ii)=dv(2)/Y(ii);
    b=(u1(ii)^2+v1(ii)^2)^(1/2);
    u1(ii)=u1(ii)/b;
    v1(ii)=v1(ii)/b;
end

[ts,xs]=ode45(dx,[0 5],[29.5;15]);
figure;
hold on
quiver(X,Y,u,v,'r')
plot(xs(:,1),xs(:,2),'g','Linewidth',4)
xlabel('mean')
ylabel('variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
figure;
quiver(X,Y,u1,v1,'b')

xlabel('mean')
ylabel('variance')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

