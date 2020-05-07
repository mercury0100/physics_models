Assignment 1 d)

% Clear memory and only show a few digits
clear all; format short;

% Dimensionalisation & parameters
C=0.1;
m=3;
rho=1.2;
radius=4.5e-2;
A=pi*radius^2;
G=9.8; % Acceleration due to gravity (m/s^2)
LS=1.0; % Choice for scaling length (m)
TS=sqrt(LS/G); % Choice for scale for time (s)

% Prompt user for time step, initial speed and angle
%tau_s=input('Enter time step in s: ');
speed_m=input('Enter initial speed in m/s: ');
%angle=input('Enter initial angle in degrees: ');
tau_s=0.00001;
%speed_m=300;
angle=0;

% Convert to radians
angle=angle*pi/180;

% Non-dimensionalise time step and initial speed
tau=tau_s/TS;
speed=speed_m/(LS/TS);

% Row vectors for non-dimensional position and velocity, and also 
% the exact position
pos=[0 20];
pose=pos;
vel=speed*[cos(angle) sin(angle)];

Del=0.5*rho*C*A*LS/m;

% Midpoint method integration
n=0;
while pos(2)>=0
  
  % Store positions for plotting
  n=n+1;
  x(n)=pos(1); y(n)=pos(2);
  xe(n)=pose(1); ye(n)=pose(2);

  % One step of the midpoint method
  accel=-Del*norm(vel)^2*vel-[0 1];
  
  vel_next=vel+tau*accel;
  pos_next=pos+0.5*tau*(vel+vel_next);
  
  pos=pos_next;
  vel=vel_next;
end

%Calculate time of flight
dydt=(pos(2)-y(n))/tau_s;
tof=tau_s*(n)+y(n)/dydt;
disp(['Time of flight (s): ',num2str(tof)]);

% Plot the midpoint method trajectory and the analytic trajectory
plot(LS*x,LS*y,'o',LS*xe,LS*ye,'-')
xlabel('Distance (m)')
ylabel('Height (m)')

% Linear interpolation to estimate the range of the projectile,
% with conversion to m
range=pos(1)-pos(2)*(pos(1)-x(end))/(pos(2)-y(end));
range_m=LS*range;
disp(['Range (m): ',num2str(range_m)])
