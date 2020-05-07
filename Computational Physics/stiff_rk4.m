% stiff_rk4.m
% Attempt to solve a stiff problem using RK4

% Clear memory and show only a few digits
clear all; format short;

% Global variable (visible inside RHS function)
global gamma;

% Stiffness parameter and initial condition
gamma=10; 
x=[1 0];

% Integration time, number of steps, and time step
T=50; nsteps=360; tau=T/nsteps;

% Eigenvalues
lam1=-gamma+sqrt(gamma^2-1);
lam2=-gamma-sqrt(gamma^2-1);

% Initial values in arrays for plotting
pos(1)=x(1);
t(1)=0;
pos_an(1)=x(1);
  
% Coefficients in analytic solution
A=(x(2)-lam2*x(1))/(lam1-lam2);
B=(lam1*x(1)-x(2))/(lam1-lam2);

% RK4 method integration
for n=1:nsteps
  
  % One step of RK4
  f1=rhs_stiff(x);
  f2=rhs_stiff(x+0.5*tau*f1);
  f3=rhs_stiff(x+0.5*tau*f2);
  f4=rhs_stiff(x+tau*f3);
  x=x+tau*(f1+2*f2+2*f3+f4)/6;
  
  % Fill in arrays with numerical results and analytic solution
  pos(n+1)=x(1);
  t(n+1)=n*tau;
  pos_an(n+1)=A*exp(lam1*t(n+1))+B*exp(lam2*t(n+1));
  
end 

% Plot numerical results and the analytic solution
plot(t,pos,'ro-',t,pos_an,'-');
axis([0 T min(pos_an) max(pos_an)]);
title(['Number of steps: ',num2str(nsteps)]); xlabel('t'); ylabel('x');
anno=legend('Numerical solution','Analytic solution');
set (anno,'Box','off')
