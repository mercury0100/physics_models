% rhs_stiff.m
% Function to evaluate the right hand side of the stiff ODE 
% considered in the lectures.
%
% Inputs:
% x - the current value of the dependent variable. For this case 
%   x = [x,v] with v=dx/dt. (The general notation becomes slightly
%   confusing for this problem.)
% Outputs:
% rhs - a row vector representing the value of the right hand side 
%   of the ODEs. Here rhs=[v,-gamma*v-x].

function rhs=rhs_stiff(x) 

% Stiffness parameter
global gamma;

rhs(1)=x(2);
rhs(2)=-2*gamma*x(2)-x(1);