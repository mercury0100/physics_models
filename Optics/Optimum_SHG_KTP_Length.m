clear all
format short

%Global variables
L=1e-3;
go=25;
R=0.998;
gamma=2;
l=[0:0.1:10]*1e-3;

P_SHG=(l.^2).*(L*go./(gamma*L+0.5*(1-exp(-2*gamma.*l)))-1).^2;

plot(l,P_SHG)
title('Output power vs Crystal length')