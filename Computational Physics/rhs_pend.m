% rhs_pend.m
% Function to evaluate the right hand side of the coupled (non-
% dimensional) ODEs describing the nonlinear pendulum
%
% Inputs:
% x - the current value of the dependent variable. For the pendulum 
%   ODEs x = [theta omega] where theta is the angle and omega is the
%   angular velocity.
% Outputs:
% rhs - a row vector representing the value of the right hand side 
%   of the ODEs. Specifically, rhs=[omega -4*pi^2*sin(theta)].

function rhs=rhs_pend(x) 

theta=x(1);
omega=x(2);
rhs(1)=omega;
rhs(2)=-4*pi^2*sin(theta);