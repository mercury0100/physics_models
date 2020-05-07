% Velocity-Verlet method to solve projectile with a magnus force

% Clear memory and only show a few digits
clear all; format short;

% Dimensionalisation parameters
G=9.8; % Acceleration due to gravity (m/s^2)
LS=1.0; % Choice for scaling length (m)
TS=sqrt(LS/G); % Choice for scale for time (s)

% Prompt user for time step, initial speed and angle
tau_s=input('Enter time step in s: ');
speed_m=input('Enter initial speed in m/s: ');
angle=input('Enter initial angle in degrees: ');
w=input('Enter angular momentum in degrees: ');
rho=1.2;
Rad=0.045;
m=0.045;



% Convert to radians
angle=angle*pi/180;
w=w*pi*TS/180;

% Non-dimensionalise time step and initial speed
tau=tau_s/TS;
speed=speed_m/(LS/TS);

Lambda=(w*rho*pi*Rad^3)/(2*m);

% Row vectors for non-dimensional position and velocity, and also 
% the exact position
pos=[0 30];
vel=speed*[cos(angle) sin(angle)];

n=0;
while pos(2)>=0
   
    n=n+1;
    
  % Store position and time for plotting: time step n
    x(n)=pos(1);
    y(n)=pos(2);
 
  % Verlet Method: position
    accel=Lambda*(vel*[0 -1;1 0])+[0 -1];
    
    acc(n)=accel(1);
    bcc(n)=accel(2);
    
    pos_next=pos+tau*vel+0.5*tau^2*accel;
    
  % Euler method: velocity  
    vel_next=vel+tau*accel;
  
  pos=pos_next;
  vel=vel_next;
  
end

figure(1)
plot(acc,bcc);

% Plot the Verlet's method trajectory
figure(2)
plot(LS*x,LS*y,'o')
xlabel('Distance (m)')
ylabel('Height (m)')

% Linear interpolation to estimate the range of the projectile,
% with conversion to m
range=pos(1)-pos(2)*(pos(1)-x(end))/(pos(2)-y(end));
range_m=LS*range;
disp(['Range (m): ',num2str(range_m)])
