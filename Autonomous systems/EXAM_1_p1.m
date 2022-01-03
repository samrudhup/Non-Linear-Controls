%Generic ODE45 code
%create time domain
tspan=[0 15];
% For loop for repeated tests
for ii = 1:100
    %Random system constants
    b=10*rand;
    c=10*rand+(1/b);
    a=10*rand+(1/b);
    %Initial conditions
    x0=[10*rand,10*rand,10*rand];
    %Use an ode solver ie ode45
    [t,y]=ode45(@(t,x)odefcn(t,x,[a,b,c]),tspan,x0);
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
end

function dxdt=odefcn(t,x,constants)
dxdt=zeros(3,1);
a=constants(1);
b=constants(2);
c=constants(3);
%ODES
dxdt(1)=-(a*x(1))+x(2);
dxdt(2)=-(b*x(2))+x(3);
dxdt(3)=-(c*x(3));
end
 