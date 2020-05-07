% numerical_error.m
% Demonstration of rounding error and truncation error in estimating 
% a derivative via [f(x+h)-f(x)]/h for f(x)=x^2 at x=1 (after Garcia 
% Section 1.5, p.28)

clear all; % Clear memory

% Logarithmically spaced values of h
NH=21;
h=logspace(-20,0,NH);

% An inline function representing the FD approximation to the 
% derivative. The following functions evaluate [f(x+h)-f(h)]/h 
% with f(x) = x^2, as well as the analytic simplification of
% this (2x+h).
fdapprox=@(x,h) ((x+h).^2-x.^2)./h;
fdapprox2=@(x,h) 2*x+h;

% The absolute error for the given h values at x = 1
error=abs(2-fdapprox(1,h));
error2=abs(2-fdapprox2(1,h));

% Is the error identically zero anywhere?
ii=find(error == 0);
NZ=length(ii);
ii2=find(error2 == 0);
NZ2=length(ii2);

% Plot the error
loglog(h,error,'o',h,error2,'ro',h,eps*ones(1,NH),'k');
axis([min(h) max(h) 1.e-1*eps 5*max(error)])
xlabel('h');
ylabel('Absolute error');

% Add annotation
handle=legend('[f(x+h)-f(x)]/h','2x+h','eps');
set (handle,'Box','off','Location','EastOutside')

disp('The error in using 2x+h is exactly zero for:');
disp(['h = ',num2str(h(ii2))]);

