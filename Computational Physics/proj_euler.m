% proj_euler.m
% Simple projectile motion, using Euler's method

% Clear memory and only show a few digits
clear all; format short;

% Dimensionalisation parameters
G=9.8; % Acceleration due to gravity (m/s^2)
LS=1.0; % Choice for scaling length (m)
TS=sqrt(LS/G); % Choice for scale for time (s)

% Non-dimensional timestep 
tau=0.1; 

% Prompt user for initial speed and angle
speed_m=input('Enter initial speed in m/s: ');
angle=input('Enter initial angle in degrees: ');

% Convert angle to radians
angle=angle*pi/180;

% Non-dimensionalise initial speed
speed=speed_m/(LS/TS);

% Row vectors for non-dimensional position and velocity
pos=[0 0];                         
vel=speed*[cos(angle) sin(angle)];

% Euler's method integration
n=0;
while pos(2)>=0
  
  % Store position for plotting
  n=n+1;
  x(n)=pos(1); y(n)=pos(2);

  % One step of Euler's method
  pos=pos+tau*vel;
  vel=vel+tau*[0 -1];
end

% Plot the trajectory as circles
plot(LS*x,LS*y,'o')
xlabel('Distance (m)')
ylabel('Height (m)')

% Linear interpolation to estimate the range of the projectile,
% with conversion to m
range=pos(1)-pos(2)*(pos(1)-x(end))/(pos(2)-y(end));
range_m=LS*range;
disp(['Range (m): ',num2str(range_m)])

% Analytic expression for range
an_range_m=speed_m^2*sin(2*angle)/G;
disp(['Analytic value for range (m): ',num2str(an_range_m)])
