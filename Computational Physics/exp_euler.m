% exp_euler.m
% Integrate dx/dt=x using Euler's method

% Clear memory and show only a few digits
clear all; format short;

% Number of steps and time step
nsteps=10; tau=1/nsteps;

% Initial values
t=0; x=1;               

% Show values of independent and dependent variables
disp(''); disp('      t       x');
fprintf('%7.5g %7.5g\n',t,x);

% Euler's method integration
for n=1:nsteps
  
  % One step of Euler
  f=rhs_exp(x);
  x=x+tau*f;
  t=t+tau;
  
  % Show values of independent and dependent variables
  fprintf('%7.5g %7.5g\n',t,x);
end 

% Display percentage error
disp(['Percentage error: ',num2str(100*abs(x-exp(1))/exp(1))]);
