% EM_rk4mod.m
% Solve the dynamical equations for a charged particle in an
% electromagnetic field

%global theta
%global alpha

% Clear memory and show only a few digits
clear all; format short;

% Name of function evaluating RHS of ODEs
rhs_string='rhs_EM';  

% Paramaters
theta = pi/4;
alpha = 1/2;
R=1;

% Initial conditions
x=[0 0 0 0 0 0];% alpha theta];

% Initial values in arrays for plotting
posx(1)=x(1);
posy(1)=x(2);
posz(1)=x(3);
t(1)=0;
pos_an(1)=x(1);
posy_an(1)=x(2);
posz_an(1)=x(3);

% Time step, total integration time and number of time steps
tau=0.05; T=10;
nsteps=ceil(T/tau);

figure(1);
% Fourth-order Runge-Kutta integration
for n=1:nsteps
    
  % Time
  time=(n-1)*tau;
  
% One step of RK4
  f1=rhs_EM(x);
  f2=rhs_EM(x+0.5*tau*f1);
  f3=rhs_EM(x+0.5*tau*f2);
  f4=rhs_EM(x+tau*f3);
  x=x+tau*(f1+2*f2+2*f3+f4)/6;
  
  % Co-ordinates of the particle
  xpart=x(1);
  ypart=x(2);
  zpart=x(3);
  

   % Fill in arrays with numerical results and analytic solution
  posx(n+1)=x(1);
  posy(n+1)=x(2);
  posz(n+1)=x(3);
  t(n+1)=n*tau;
  posx_an(n+1)=R*sin(time);
  posy_an(n+1)=R*cos(time)-alpha*sin(theta*time);
  posz_an(n+1)=0.5*alpha*cos(theta*time^2);
  
  % Animate the pendulum motion
  plot3(xpart,ypart,zpart,'o-')
  title(['Time: ',num2str(t(n)+tau,'%5.3f')])
  axis equal
  axis([-1 1 -1 1])
  pause(0.05);
  drawnow;
  
  
end 

% Plot numerical results and the analytic solution
plot3(posx,posy,posz,'ro-',posx_an,posy_an,posz_an,'-');
%axis([0 T min(pos_an) max(pos_an)]);
title(['Number of steps: ',num2str(nsteps)]); xlabel('t'); ylabel('x');
anno=legend('Numerical solution','Analytic solution');
set (anno,'Box','off')
