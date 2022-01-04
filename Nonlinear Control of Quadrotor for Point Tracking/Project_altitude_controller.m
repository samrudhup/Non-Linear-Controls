clc; clear;
close all;

%defining the constants and their conditions
g= 9.81;          %gravity
zd=5*rand + 0.01;
tspan = [0 15]; 
m=5*rand + 0.01;%mass
k7=5*rand + 0.01;%control parameter
k8=5*rand + 0.01;%control parameter
phi=5*rand + 0.01;
theta=5*rand + 0.01;


for ii=1:100
    theta=5*rand + 0.01;
    phi=5*rand + 0.01;
    e7_0=5*rand + 0.01;
    e8_0=5*rand + 0.01;
    z_dot0=0;
    V0=0;
    y0=[e7_0,e8_0,z_dot0,V0];
    [t,y]   = ode45(@(t,y) alt(t,y,m,k7,k8,g,zd,phi,theta),tspan, y0);
    figure(1)
    plot(t,y(:,1));
    title('Resulting Trajectories of the 1st State w. 100 Simulations')
    xlabel('Time(s)')
    ylabel('Error, e7(t)')
    hold on
    figure(2)
    plot(t,y(:,2));
    title('Resulting Trajectories of the 2nd State w. 100 Simulations')
    xlabel('Time(s)')
    ylabel('Error, e8(t)')
    hold on
    figure(3)
    plot(t,y(:,3));
    title('Resulting Trajectories of the 3rd State w. 100 Simulations')
    xlabel('Time(s)')
    ylabel('z_ddot(t)')
    hold on
    
    
    figure(4)
    plot(t,y(:,4));
    title('Resulting Trajectories of the 5th State w. 100 Simulations')
    xlabel('Time(s)')
    ylabel('Lyapunov function')
    hold on
end
function dydt = alt(t,y,m,k7,k8,g,zd,phi,theta)
    dydt=zeros(4,1);%create an empty matrix 
    e7= y(1);    %tracking error
    e8= y(2);
    z_dot= y(3);     %estimated altitude of the quadrotor
    V=y(4);
    zd_ddot=0;  %desired velocity of the altitude is assumed to be zero
    
    z=e7+zd;
    %Design U1
    U1=(m*g + ((k7^2)-1)*m*e7 - (k7+k8)*e8*m)/(cos(phi)*cos(theta));
    
    e7_dot= e8 - k7*e7;
    e8_dot=(cos(phi)*cos(theta)*(U1/m)-g-zd_ddot+k7*(e8-k7*e7));
    z_ddot=dydt(2)-k7*e7_dot+zd_ddot;
    
    V_dot=e7*e7_dot+e8*e8_dot;
    
    
    dydt(1)=e7_dot;
    dydt(2)=e8_dot;
    dydt(3)=z_ddot;
    dydt(4)=V_dot;
end
    
    