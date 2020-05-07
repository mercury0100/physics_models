% lab02_q2b.m
% Solution file for Lab 2 Q2 (b)

% Clear memory and only show a few digits
clear all; format short;

% Dimensionalisation parameters
G=9.8; % Acceleration due to gravity (m/s^2)
LS=1.0; % Choice for scaling length (m)
TS=sqrt(LS/G); % Choice for scale for time (s)

% Prompt user for time step, initial speed and angle
%tau_s=input('Enter time step in s: ');
%speed_m=input('Enter initial speed in m/s: ');
%angle=input('Enter initial angle in degrees: ');
tau_s=0.01;
speed_m=10;
angle=60;

% Convert to radians
angle=angle*pi/180;

% Non-dimensionalise time step and initial speed
tau=tau_s/TS;
speed=speed_m/(LS/TS);

% Row vectors for non-dimensional position and velocity, and also 
% the exact position
pos=[0 0];
pose=pos;
vel=speed*[cos(angle) sin(angle)];
vel1=vel;

% Midpoint method integration
n=0;
while pos(2)>=0
  
  % Store positions for plotting
  n=n+1;
  x(n)=pos(1); y(n)=pos(2);
  xe(n)=pose(1); ye(n)=pose(2);

  % One step of the midpoint method
  velp=vel; % Store value of velocity
  vel=velp+tau*[0 -1];
  pos=pos+0.5*tau*(vel+velp);
  
  % Analytic solution
  time=n*tau;
  pose=vel1*time+0.5*[0 -1]*time^2;
end

% Plot the midpoint method trajectory and the analytic trajectory
plot(LS*x,LS*y,'o',LS*xe,LS*ye,'-')
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

% Percentage error
err=(range_m/an_range_m-1)*100;
disp(['Percentage error in range: ',num2str(err)])
