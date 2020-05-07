% kepler_midpoint.m
% Motion under a central force, using the midpoint method

% Clear memory and only show a few digits
clear all; format short;

% Time step (non-dim.)
tau=0.05;

% Initial position (non-dim.) - this should be fixed
pos=[1 0];

% Initial velocity (non-dim.) - vary the y-component
vel=[0 1];

% Total integration time
T=4*pi;

% Number of integration steps
steps=ceil(T/tau);

% Plot only 100 frames in total
skip=ceil(steps/100);

% Calculate trajectory from analytic solution. See Appendix
% of Week 2 lecture 1.
[xan,yan]=kepler_analytic(vel,T);

% Midpoint method integration
figure(1);
for n=1:steps

  % Store position and time for plotting: time step n
  x(n)=pos(1); y(n)=pos(2); time(n)=(n-1)*tau;
  
  % Plot numerical and analytic solution
  if rem(n,skip)==0
    plot(x,y,'g-',pos(1),pos(2),'ko',xan,yan,'b',0,0,'ro')
    title(['Time: ',num2str(time(n))]);
    xlabel('x');
    ylabel('y');
    axis equal; % Preserve aspect ratio
    drawnow; % Draw immediately
  end

  % Calculate radial position, speed and acceleration at step n
  r=norm(pos); speed=norm(vel); accel=-pos/r^3;
 
  % Calculate total energy at step n and store
  energy(n)=0.5*speed^2-1/r;
  
  % One step of the midpoint method
  velp=vel; % Store old velocity
  vel=vel+tau*accel;
  pos=pos+0.5*tau*(velp+vel);
  
end

% Plot energy versus time
figure(2);
plot(time,energy);
xlabel('Time (non-dim.)');
ylabel('Total energy (non-dim.)');
