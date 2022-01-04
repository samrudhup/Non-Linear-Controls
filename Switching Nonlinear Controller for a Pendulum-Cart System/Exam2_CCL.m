clc;clear;
%defining all the constants
g=9.81;%gravity
a_hat_x = rand + 0.5; %estimated input coefficient
a_hat_phi=rand + 0.5; 
m_hat_x = 0.4; %estimated mass
m_hat_phi = 0.5;
c_hat_x = 0.6; %estimated damper coefficient
c_hat_phi=  0.6;
%k_hat = 5*rand + 0.01; %estimated spring coefficient
%phid   = -5 + (5+5)*rand; %tilt angle has no restriction, could be negative
alpha = 1.4; %design variable, (must be greater than zero)
beta  = 1.8; %design variable, (must be greater than zero)
%xd    = -1 + (1+1)*rand;%desired location of the mass, positive constant
%theta_hat_x = [m_hat_x; m_hat_phi; m_hat_phi*l; c_hat_x];
%theta_hat_phi=[m_hat_phi*l; m_hat_phi*l^2; c_hat_phi];
lambda = rand + 0.5; %design variable, (must be greater than zero)
gamma_a = rand + 0.5;
tspan = (0:0.01:30); 
freq_phi=0.75;
freq_x = 0.7;
l=0.8;
xd_bar    = 0.8;
phid_bar  = pi;
V_0=0;
    
for ii=1:100
    e_x_0=rand + 0.5;
    e_phi_0=rand + 0.5;
    r_x_0=rand + 0.5;
    r_phi_0=rand + 0.5;
    a_til_x_0=rand + 0.5;
    a_til_phi_0=rand + 0.5;
    m_til_x_0=0.5+rand;
    m_til_phi_0=0.5+rand;
    c_til_x_0=0.5+rand;
    c_til_phi_0=rand+0.5;
    u_til_x_0=0;
    u_til_phi_0=0;
    x_dot_0=0;
    phi_dot_0=0;
   
    y0=[r_x_0, u_til_x_0, e_x_0, m_til_x_0, c_til_x_0, a_til_x_0, x_dot_0, phi_dot_0, r_phi_0, u_til_phi_0, e_phi_0, m_til_phi_0, c_til_phi_0, a_til_phi_0, V_0];
    
    [t,y]   = ode45(@(t,y) ode45ccl(t,y,a_hat_x, m_hat_x, c_hat_x, g, alpha, beta, lambda, gamma_a, freq_x, a_hat_phi, m_hat_phi, c_hat_phi, freq_phi, l, xd_bar, phid_bar), tspan, y0);
   
    %linear plot for each state
    figure(1)
    plot(t,(y(:,1)));
    title('Filter tracking error r(t) of x')
    xlabel('Time(s)')
    ylabel('Filered tracking error r(t) of x')
	hold on
    
    figure(2)
    plot(t,y(:,2),t,y(:,6),t,y(:,10),t,y(:,14));
    title('Norms of the tracking errors u and a')
    xlabel('Time(s)')
    ylabel('vector norm')
    hold on
    
    figure(3)
    plot(t,y(:,3));
    title('Resulting Trajectories of the 3rd State w. 100 Simulations')
    xlabel('Time(s)')
    ylabel('Tracking Error, e(t) of x')
    hold on
    
    figure(4)
    subplot(2,1,1)
    plot(t,y(:,4));
    title('Resulting Trajectories of mass of x')
    hold on
    subplot(2,1,2)
    plot(t,y(:,5));
    title('Resulting Trajectories of mass of x')
    hold on
    xlabel('Time(s)')
    ylabel('Parametric errors of the collar')
    hold on
    
    figure(5)
    plot(t,y(:,6));
    title('Resulting Trajectories of the 5th State w. 100 Simulations')
    xlabel('Time(s)')
    ylabel('Parametric Input Coefficient Error, a_{til}(t) of x')
    hold on
    
    figure(6)
    plot(t,y(:,9));
    title('Resulting Trajectories of the 1st State w. 100 Simulations')
    xlabel('Time(s)')
    ylabel('Filtered Tracking Error, r(t) of phi')
    hold on
    
    figure(7)
    plot(t,y(:,10));
    title('Resulting Trajectories of the 2nd State w. 100 Simulations')
    xlabel('Time(s)')
    ylabel('Backstepping Error, u_{til}(t) of phi')
    hold on
    
    figure(8)
    plot(t,y(:,11));
    title('Resulting Trajectories of the 3rd State w. 100 Simulations')
    xlabel('Time(s)')
    ylabel('Tracking Error, e(t) of phi')
    hold on
    
    figure(9)
    subplot(2,1,1)
    plot(t,y(:,12));
    title('Resulting Trajectories of mass for phi')
    hold on
    subplot(2,1,2)
    plot(t,y(:,13));
    title('Resulting Trajectories of mass for phi')
    hold on
    xlabel('Time(s)')
    ylabel('Parametric errors of the pendulum')
    hold on
    
    figure(10)
    plot(t,y(:,14));
    title('Resulting Trajectories of the 5th State w. 100 Simulations')
    xlabel('Time(s)')
    ylabel('Parametric Input Coefficient Error, a_{til}(t) of phi')
    hold on
end

function dydt=ode45ccl(t,y,a_hat_x, m_hat_x, c_hat_x, g, alpha, beta, lambda, gamma_a, freq_x,a_hat_phi, m_hat_phi, c_hat_phi, freq_phi, l,xd_bar, phid_bar)
dydt=zeros(15,1);
r_x     =y(1);
u_til_x =y(2);
e_x     =y(3);
m_til_x =y(4);
c_til_x =y(5);
a_til_x =y(6);
x_dot   =y(7);
phi_dot =y(8);
r_phi   =y(9);
u_til_phi=y(10);
e_phi   =y(11);
m_til_phi=y(12);
c_til_phi=y(13);
a_til_phi=y(14);



%DESIGNS:(Ud & mu)


xd        = xd_bar*sin(2*pi*freq_x*t);
xd_dot    = 2*pi*freq_x*xd_bar*cos(2*pi*freq_x*t);
xd_ddot   = -((2*pi*freq_x)^2)*xd_bar*sin(2*pi*freq_x*t);
xd_dddot  = -((2*pi*freq_x)^3)*xd_bar*cos(2*pi*freq_x*t);

phid      = phid_bar*sin(2*pi*freq_phi*t);
phid_dot  = 2*pi*freq_phi*phid_bar*cos(2*pi*freq_phi*t);
phid_ddot   = -((2*pi*freq_phi)^2)*phid_bar*sin(2*pi*freq_phi*t);
phid_dddot  = -((2*pi*freq_phi)^3)*phid_bar*cos(2*pi*freq_phi*t);

x         = xd - e_x;   %position of the mass
e_dot_x   = -x_dot;
e_ddot_x  = dydt(1)-alpha*dydt(3);
x_ddot    = -e_ddot_x;

phi       = phid - e_phi; % angular position of the arm
e_dot_phi = -phi_dot;
e_ddot_phi= dydt(9)-alpha*dydt(11);
phi_ddot  = -e_ddot_phi;

theta_hat_x = [m_hat_x; m_hat_phi; m_hat_phi*l; c_hat_x];
theta_hat_phi=[m_hat_phi*l; m_hat_phi*(l^2); c_hat_phi];

%designed ud for x
Y_x     = [ xd_ddot+alpha*dydt(3); (alpha*dydt(3)*(sin(phi)^2)) + (sin(phi)^2)*xd_ddot - g*cos(phi)*sin(phi) + (sin(phi)*cos(phi)*phi_dot*r_x); (-(phi_dot)^2)*sin(phi); x_dot];
Y_dot_x = [ xd_dddot+alpha*e_ddot_x; (2*sin(phi)*cos(phi)*phi_dot*xd_ddot) + ((sin(phi)^2)*xd_dddot) + g*phi_dot*(sin(phi)^2) - g*phi_dot*(cos(phi)^2) + r_x*sin(phi)*cos(phi)*dydt(14) + r_x*(cos(phi)^2)*(phi_dot^2) - (sin(phi)^2)*r_x*(phi_dot^2) + dydt(1)*sin(phi)*cos(phi)*phi_dot + alpha*dydt(3)*2*sin(phi)*cos(phi)*phi_dot + alpha*(sin(phi)^2)*e_ddot_x; -cos(phi)*(phi_dot^3) - 2*sin(phi)*phi_dot*dydt(14); x_ddot];

ud_x    = (Y_x')*theta_hat_x + beta*r_x + e_x;
u_x = ud_x-u_til_x;
m_x = m_hat_x + m_til_x; 
c_x = c_hat_x + c_til_x; 
a_x = a_hat_x + a_til_x;

ud_dot_x= (Y_dot_x.')*theta_hat_x + beta*dydt(1) + dydt(3);
mu_x    = a_hat_x*u_x + ud_dot_x + r_x + lambda*u_til_x;
x_dddot   = (-((c_hat_x + c_til_x)*x_ddot)-(a_hat_x + a_til_x)*(ud_x-u_til_x)+mu_x)/((m_hat_x + m_til_x)+(m_hat_phi + m_til_phi)*(sin(phi)^2));

%designed ud for phi
Y_phi   = [ g*sin(phi) + cos(phi)*x_ddot; phid_ddot + alpha*e_dot_phi; phi_dot];
Y_dot_phi=[(g*cos(phi)*phi_dot + cos(phi)*x_dddot - sin(phi)*phi_dot*dydt(7)); phid_dddot+alpha*e_ddot_phi; dydt(8)];

ud_phi  = (Y_phi.')*theta_hat_phi + beta*r_phi + e_phi;
ud_dot_phi= (Y_dot_phi.')*theta_hat_phi + beta*dydt(9) + dydt(11);


u_phi = ud_phi-u_til_phi;
m_phi = m_hat_phi + m_til_phi; 
c_phi = c_hat_phi + c_til_phi; 
a_phi = a_hat_phi + a_til_phi;


mu_phi  = a_hat_phi*u_phi+ud_dot_phi+r_phi+lambda*u_til_phi;

Gamma1=eye(4);
Gamma2=eye(3);

%defining outputs
%for x
r_dot_x =(m_x*xd_ddot + m_phi*(sin(phi)^2)*xd_ddot + (c_x*x_dot) - m_phi*g*cos(phi)*sin(phi) - m_phi*l*(phi_dot^2)*sin(phi) + u_til_x - ud_x + m_x*alpha*dydt(3) + m_phi*(sin(phi)^2)*alpha*e_dot_x)/(m_x + m_phi*(sin(phi)^2));
u_til_dot_x=(ud_dot_x + a_x*(ud_x-u_til_x)-mu_x);
e_dot_x  = r_x-alpha*e_x;
theta_til_dot_x = Gamma1 .* Y_x' .* r_x;
m_til_dot_x = theta_til_dot_x(1);
c_til_dot_x = theta_til_dot_x(4);
m_til_dot_phi = theta_til_dot_x(2);
a_til_dot_x=- gamma_a*(ud_x*u_til_x - u_til_x^2);
%for phi
r_dot_phi =((m_phi*(l^2)*phid_ddot)+ c_phi*phi_dot + m_phi*l*g*sin(phi)+ m_phi*l*cos(phi)*dydt(7)+ u_til_phi - ud_phi + alpha*m_phi*(l^2)*dydt(11))/(m_phi*(l^2));
u_til_dot_phi =(ud_dot_phi + a_phi*u_phi - mu_phi);
e_dot_phi =r_phi - alpha*e_phi;
theta_til_dot_phi = Gamma2 .* Y_phi' .* r_phi;
c_til_dot_phi = theta_til_dot_phi(3);
a_til_dot_phi=- gamma_a*(ud_phi*u_til_phi - (u_til_phi^2));

V_dot_x=e_x*dydt(3)+m_x + m_phi*(sin(phi)^2)*r_x*dydt(1)+(m_phi*cos(phi)*sin(phi)*phi_dot*r_x^2)+u_til_x*dydt(10)+(a_til_x*dydt(6))/gamma_a;

dydt      = [r_dot_x; u_til_dot_x; e_dot_x; m_til_dot_x; c_til_dot_x; a_til_dot_x; x_ddot; phi_ddot; r_dot_phi; u_til_dot_phi; e_dot_phi; m_til_dot_phi; c_til_dot_phi; a_til_dot_phi; V_dot_x]; 
end
