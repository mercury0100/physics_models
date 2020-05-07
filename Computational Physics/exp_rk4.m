% exp_rk4.m
% Integrate dx/dt=x using RK4

% Clear memory and show only a few digits
clear all; format short;

% Number of steps and time step
nsteps=10; tau=1/nsteps;

% Initial values
t=0; x=1;               

% Show values of independent and dependent variables
disp(''); disp('      t       x');
fprintf('%7.5g %7.5g\n',t,x);

% RK4 method integration
for n=1:nsteps
  
  % One step of RK4
  f1=rhs_exp(x);
  f2=rhs_exp(x+0.5*tau*f1);
  f3=rhs_exp(x+0.5*tau*f2);
  f4=rhs_exp(x+tau*f3);
  x=x+tau*(f1+2*f2+2*f3+f4)/6;
  t=t+tau;
  
  % Show values of independent and dependent variables
  fprintf('%7.5g %7.5g\n',t,x);
end 

% Display percentage error
disp(['Percentage error: ',num2str(100*abs(x-exp(1))/exp(1))]);
