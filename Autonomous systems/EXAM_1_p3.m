%Generic ODE45 code
%create time domain
clc
clear
tspan=[0 20];
% For loop for repeated tests
for ii = 1:100
    %Random system constants
    a=rand+0.5;
    c=rand+0.5;
    m=rand+0.5;
    g=9.81;
    k=rand+0.5;
    alpha=rand+0.5;
    beta=rand+0.5;
    phi=rand+0.5;
    mu=0;
    xd=rand;
    %Initial conditions
    x0=[rand+0.2,rand+0.2,rand+0.2,rand+0.2];
    %Use an ode solver ie ode45
    [t,y]=ode45(@(t,x)odefcn(t,x,a,c,m,g,k,alpha,beta,phi,xd),tspan,x0);
    %Plot results
    figure(1)
    plot(t,y(:,1))
    xlabel('t');
    ylabel('x1_dot');
    hold on
    figure(2)
    plot(t,y(:,2))
    xlabel('t');
    ylabel('x2_dot');
    hold on
    figure(3)
    plot(t,y(:,3))
    xlabel('t');
    ylabel('x3_dot');
    hold on
    figure(4)
    plot(t,y(:,4))
    xlabel('t');
    ylabel('x4_dot');
    hold on
    
end

function dxdt=odefcn(t,x,a,c,m,g,k,alpha,beta,phi,xd)
dxdt=zeros(4,1);

%ODES
%Ud_dot=0;
% Ud=0
Ud=-c*dxdt(1)+k*(xd-x(1))-m*g*sin(phi)+x(3)+m*alpha*dxdt(1)+beta*x(2)+x(1);
Ud_dot=dxdt(2)*(-c+m*alpha+beta)+dxdt(1)*(c*alpha-k-m*alpha^2+1)+dxdt(3);

mu=Ud_dot+a+Ud;

dxdt(1)=x(2)-alpha*x(1);

dxdt(2)=(-c*dxdt(1)+k*xd-k*x(1)-(m*g*sin(phi))+x(3)-Ud)/m+alpha*dxdt(1);

dxdt(3)=Ud_dot+a*Ud-a*x(3)-mu;


dxdt(4)=mu;
end
 