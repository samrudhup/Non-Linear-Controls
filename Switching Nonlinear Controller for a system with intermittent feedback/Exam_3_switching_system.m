clc; clear;
%Defining all the constants
m_hat=rand+0.5; %estimated mass 
c_hat=rand+0.5; %estimated friction coefficient
alpha=rand+0.5;
beta_r=rand+0.5;
beta_d=rand+0.5;
lambda = eye(2); %design variable, (must be greater than zero)
gamma_a = rand + 0.5;
tspan = (0:0.01:40);
freq=0.003;
R=100;
D=normrnd(0,3,2,1);
%new_dot_m is only available when every max dwell time units

for ii= 1:1
    r_hat_x_0=rand+0.5;
    r_hat_y_0=rand+0.5;
    m_til_0=rand+0.5;
    c_til_0=rand+0.5;
    tau_til_0=0;
    new_til_0=[rand+0.5 rand+0.5].';
    new_0=[ 105 -5].';
    x_dot_0=0;
    y_dot_0=0;
    new_dot_0=[ x_dot_0 y_dot_0].';
    new_ddot_0=[ 0 0].';
    e_hat_x_0=rand+0.5;
    e_hat_y_0=rand+0.5;
    x_axis=0;
    y_axis=0;
    
    
     y0=[r_hat_x_0, r_hat_y_0, e_hat_x_0, e_hat_y_0, m_til_0, c_til_0, x_dot_0, y_dot_0];
     
     [t,y] = ode45(@(t,y) ode45cclswitch(t,y, m_hat, c_hat, alpha, beta_d, beta_r, lambda, D, R, freq), tspan, y0);
     
    %linear plot of states
    figure(1)
    plot(t,(y(:,1)));
    title('Filter tracking error r(t) of x')
    xlabel('Time(s)')
    ylabel('Filered tracking error r(t) of x')
	hold on
    
    figure(2)
    plot(t,y(:,2));
    title('Filter tracking error r(t) of y')
    xlabel('Time(s)')
    ylabel('Filter tracking error r(t) of y')
    hold on
    
    figure(3)
    plot(t,y(:,3));
    title('Tracking error e(t) of x')
    xlabel('Time(s)')
    ylabel('Tracking Error, e(t) of x')
    hold on
    
    figure(4)
    plot(t,y(:,4));
    title('Tracking Error, e(t) of y')
    xlabel('Time(s)')
    ylabel('Tracking Error, e(t) of y')
    hold on
  
    figure(5)
    plot(t,y(:,7));
    title('Resulting Trajectories of the 2nd State w. 100 Simulations')
    xlabel('Time(s)')
    ylabel('Position tracking error of x')
    hold on
    
    figure(6)
    plot(t,y(:,8));
    title('Resulting Trajectories of the 2nd State w. 100 Simulations')
    xlabel('Time(s)')
    ylabel('Position tracking error of y')
    hold on
    
    figure(7)
    plot(t,y(:,1),t,y(:,2),t, y(:,3),t,y(:,4));
    title('Resulting Trajectories of the estimation and tracking errors')
    xlabel('Time(s)')
    ylabel('Norm of errors')
    hold on
    
    figure(8)
    plot(t,y(:,5));
    title('Resulting Trajectories of the 2nd State w. 100 Simulations')
    xlabel('Time(s)')
    ylabel('Parametric error m_til')
    hold on
    
    figure(9)
    plot(t,y(:,6));
    title('Resulting Trajectories of the 2nd State w. 100 Simulations')
    xlabel('Time(s)')
    ylabel('Parametric error c_til')
    hold on
end

function dydt=ode45cclswitch(t,y, m_hat, c_hat, alpha, beta_d, beta_r, lambda, D, R, freq)
dydt=zeros(8,1);
r_hat_x=y(1);
r_hat_y=y(2);
%tau_til=y(3);
e_hat_x=y(3);
e_hat_y=y(4);
m_til=y(5);
c_til=y(6);
x_dot=y(7);
y_dot=y(8);
new_dot=[ x_dot, y_dot].';
r_hat=[r_hat_x, r_hat_y].';
e_hat=[e_hat_x, e_hat_y].';
theta_til=[m_til, c_til].';

%design tau
newd=R.*[cosd(2*pi*freq*t), sind(2*pi*freq*t)].';
newd_dot= 2*pi*freq*R.*[-sind(2*freq*pi*t), cosd(2*pi*freq*t)].';
newd_ddot= ((2*pi*freq)^2)*R.*[-cosd(2*pi*freq*t), -sind(2*pi*freq*t)].';

w_new=normrnd(0,3,2,1);
w_new_dot=normrnd(0,3,2,1);
w_new_ddot=normrnd(0,3,2,1);

new=newd-e_hat-w_new;
new_ddot=[dydt(7), dydt(8)];

num=t/(9.21);
new_0=[ 105 -5].';
if isreal(num) && rem(num,1)==0
    new_m=new_0;
else
    new_m=new+w_new;
end

new_dot_m=new_dot+w_new_dot;
new_ddot_m=new_ddot+w_new_ddot;


e_hat_dot=r_hat- alpha*(newd-new_m); %define new_dot_m

e_hat_ddot=newd_ddot - new_ddot - w_new_ddot;

new_hat=new_m;
new_hat_dot=new_dot_m;
new_hat_ddot=new_ddot_m;

theta_hat=[m_hat, c_hat].';

new_til=new-new_hat;
% x_axis=[0;x_axis;new_til(1)]
% y_axis=[0;y_axis;new_til(2)]

new_til_dot=new_dot-new_hat_dot;
Y=[newd_ddot-w_new_ddot+(alpha.*e_hat_dot), new_dot];
tau=Y*theta_hat + beta_d.*sign(r_hat) +beta_r.*r_hat+e_hat;

%define theta_til_dot
theta_til_dot=-lambda.*(Y.').*r_hat;
m_til_dot=theta_til_dot(1);
c_til_dot=theta_til_dot(2);

m=m_til+m_hat;
c=c_til+c_hat;
%defining outputs
theta=theta_til+theta_hat;
r_hat_dot=(Y*theta-tau-D)/m;
new_ddot=(-c.*new_dot+tau+D)./m;
e_hat_dot=e_hat_dot.';
dydt=[r_hat_dot(1), r_hat_dot(2), e_hat_dot(1), e_hat_dot(2), m_til_dot, c_til_dot, new_ddot(1), new_ddot(2)].';
end









    
    