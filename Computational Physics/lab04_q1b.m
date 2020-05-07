% lab04_q1b.m
% Solution code for Lab 4 Q1 (b)

% Clear memory and show only a few digits
clear all; format short;

% Prompt user for initial angle, convert to radians
theta_deg=input('Initial angle (deg): ');
theta=theta_deg*pi/180;

% Row vector with initial angle and position
x=[theta 0];

% Prompt user for time step
tau=input('Time step: ');

% Integration time and number of time steps
T=3.5;
nsteps=ceil(T/tau);

% Store arrays with time and angle
time(1)=0;
thvals(1)=x(1);

figure(1);

% Velocity-Verlet integration
for n=1:nsteps
  
  % One V-V step (in the notation of V-V, to avoid confusion)
  pos=x(1);
  vel=x(2);
  accel=-4*pi^2*sin(pos);
  pos_next=pos+tau*vel+0.5*tau^2*accel;
  accel_next=-4*pi^2*sin(pos_next);
  vel_next=vel+0.5*tau*(accel+accel_next);

  % Update vector x
  x(1)=pos_next;
  x(2)=vel_next;

  % Record time, angle in vectors
  time(n+1)=n*tau;
  thvals(n+1)=x(1);
  
  % Co-ordinates of the pendulum bar
  xpend=[0  sin(x(1))];
  ypend=[0 -cos(x(1))];
  
  % Animate the pendulum motion
  plot(xpend,ypend,'o-')
  title(['Time: ',num2str(time(n+1),'%5.3f')])
  axis equal
  axis([-1 1 -1 1])
  drawnow;
end 

% Plot angle versus time (in degrees)
thvals=180*thvals/pi;
figure(2);
plot(time,thvals,'o');
xlabel('Time (non-dim.)');
ylabel('Angle (deg.)');

