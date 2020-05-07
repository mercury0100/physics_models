% pend_rk4.m
% Solve the nonlinear pendulum problem using RK4

% Clear memory and show only a few digits
clear all; format short;

% Initial conditions
theta1=30; % Initial angle in degrees
theta=theta1*pi/180;
x=[theta 0];
vel=0

% Time step, total integration time and number of time steps
tau=0.025; T=2;
nsteps=ceil(T/tau);

figure(1);
% Fourth-order Runge-Kutta integration
for n=1:nsteps
    
  % Time
  time(n)=(n-1)*tau;

  %Velocity Verlet
  accel=-4*pi^2*sin(theta);
  theta_next=theta+tau*vel+0.5*tau^2*accel;
  accel_next=-4*pi^2*sin(theta_next);
  vel_next=vel+0.5*tau*(accel+accel_next);
  
  
  vel=vel_next;
  accel=accel_next;
  theta=theta_next;
 
  % Co-ordinates of the pendulum bar
  xpend=[0  sin(x(1))];
  ypend=[0 -cos(x(1))];
  thetas(n)=x(1);
  
  % Animate the pendulum motion
  plot(xpend,ypend,'o-')
  title(['Time: ',num2str(time(n)+tau,'%5.3f')])
  axis equal
  axis([-1 1 -1 1])
  pause(0.05);
  drawnow;
end 
