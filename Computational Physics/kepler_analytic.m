function [ xan, yan ] = kepler_analytic( vel, T )
% Calculate the analytic trajectory for the Kepler central force
% problem. 

% Calculate trajectory from analytic solution
e=norm(vel)^2-1;
a=1/(1-e);
if (e < 1)
    theta=linspace(0,2*pi,50); % Equally spaced values from 0 to 2*pi
    b=a*sqrt(1-e^2);
    xan=-a*e+a*cos(theta);
    yan=b*sin(theta);
else
    b=a*sqrt(e^2-1);
    theta_max=asinh(norm(vel)*T/b); % Limit for range of theta
    theta=linspace(-theta_max,theta_max,50);
    xan=-a*e+a*cosh(theta);
    yan=b*sinh(theta);
end

end
