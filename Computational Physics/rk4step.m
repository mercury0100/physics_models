function xout=rk4step(x,t,tau,rhs_string)

% Perform one step of RK4 for the ODE dx/dt=f(x,t). The RHS of the ODE
% is described by the function matching the supplied string rhs_string. 
%
% Routine adapted from Garcia, Numerical Methods in Physics.
%
% Inputs:
% x - current values of dependent variables
% t - current value of independent variable
% tau - step size
% rhs_string - name of function specifying the RHS of the ODE
%
% Outputs:
% xout - new value of x after a step of size tau

htau=0.5*tau;
f1=feval(rhs_string,x,t);  
th=t+htau;
xtemp=x+htau*f1;
f2=feval(rhs_string,xtemp,th);  
xtemp=x+htau*f2;
f3=feval(rhs_string,xtemp,th);
tf=t+tau;
xtemp=x+tau*f3;
f4=feval(rhs_string,xtemp,tf);
xout=x+tau*(f1+2*f2+2*f3+f4)/6;

return;