% proj_drag_euler.m
% Solve the dynamics ODEs for a cricket ball subject to gravity
% and a quadratic drag force

% Clear memory and only show a few digits
clear all; format short;

% Drag parameters
CD=0.35;           % Dimensionless drag coefficient
RHO=1.2;           % Density of air (kg/m^3)
M=0.145;           % Mass of ball (kg)
R=0.037;           % Ball radius (m)
A=pi*R^2;          % Cross-sectional area (m^2)

% Dimensionalisation parameters
G=9.8; % Acceleration due to gravity (m/s^2)
LS=1.0; % Choice for scaling length (m)
TS=sqrt(LS/G); % Choice for scale for time (s)

% Calculate non-dimensional drag parameter used in ODEs
D=0.5*CD*RHO*A*LS/M;

% Non-dimensional timestep 
tau=0.05; 

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

  % The new acceleration term
  acc=-D*norm(vel)*vel-[0 1];

  % One step of Euler's method
  pos=pos+tau*vel;
  vel=vel+tau*acc;
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
