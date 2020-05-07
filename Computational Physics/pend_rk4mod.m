% pend_rk4mod.m
% Solve the dynamic of a charged particle in an electromagnetic field using RK4, with a modular
% implementation of an RK4 step

% Clear memory and show only a few digits
clear all; format short;

% Name of function evaluating RHS of ODEs
rhs_string='rhs_pend';  

% Initial conditions
theta1=30; % Initial angle in degrees
theta=theta1*pi/180;
x=[theta 0];

% Time step, total integration time and number of time steps
tau=0.025; T=10;
nsteps=ceil(T/tau);

figure(1);
% Fourth-order Runge-Kutta integration
for n=1:nsteps
    
  % Time
  time(n)=(n-1)*tau;
  
  % Runge-Kutta step
  x=rk4step(x,time(n),tau,rhs_string);

  % Co-ordinates of the pendulum bar
  xpend=[0  sin(x(1))];
  ypend=[0 -cos(x(1))];
 
  % Animate the pendulum motion
  plot(xpend,ypend,'o-')
  title(['Time: ',num2str(time(n)+tau,'%5.3f')])
  axis equal
  axis([-1 1 -1 1])
  pause(0.05);
  drawnow;

end 
