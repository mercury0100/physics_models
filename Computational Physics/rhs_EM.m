% rhs_pend.m
% Function to evaluate the right hand side of the coupled (non-
% dimensional) ODEs describing the nonlinear pendulum
%
% Inputs:
% x - the current value of the dependent variable. For the pendulum 
%   ODEs x = [theta omega] where theta is the angle and omega is the
%   angular velocity.
% t - the current value of the independent variable.
% Outputs:
% rhs - a row vector representing the value of the right hand side 
%   of the ODEs. Specifically, rhs=[omega -4*pi^2*sin(theta)].

function rhs=rhs_EM(x) 

%global theta;
%global alpha;
alpha=1/2;
theta=0;
rhs(1)=x(4);
rhs(2)=x(5);
rhs(3)=x(6);
rhs(4)=alpha*sin(theta)+x(5);
rhs(5)=-x(4);
rhs(6)=-alpha*cos(theta);
%rhs(7)=0;
%rhs(8)=0;
