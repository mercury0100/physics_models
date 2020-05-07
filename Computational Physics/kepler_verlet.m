% kepler_verlet.m
% Motion under a central force, using the Verlet method

% Clear memory and only show a few digits
clear all; format short;

% Time step (non-dim.)
tau=0.05;

% Initial position (non-dim.) - this should be fixed
pos=[1 0];

% Initial velocity (non-dim.) - vary the y-component
vel=[0 1.3];

% Total integration time
T=1000;

% Number of integration steps
steps=ceil(T/tau);

% Plot only 100 frames in total
skip=ceil(steps/100);
                                 
% Calculate trajectory from analytic solution. See Appendix
% of Week 2 lecture 1.
[xan,yan]=kepler_analytic(vel,T);

% Verlet method integration
figure(1);
for n=1:steps
    
  % Store position and time for plotting: time step n
  x(n)=pos(1); y(n)=pos(2); time(n)=(n-1)*tau;

  % Plot numerical and analytic solution
  if rem(n,skip)==0
    plot(x,y,'g-',pos(1),pos(2),'ko',xan,yan,'b',0,0,'ro')
    title(['Time: ',num2str(time(n))]); xlabel('x'); ylabel('y');
    axis equal; % Preserve aspect ratio
    drawnow; % Draw immediately
  end

  % Calculate radial position and acceleration at step n
  r(n)=norm(pos);
  accel=-pos/r(n)^3;
 
  % Verlet Method: position - prev is step n-1, next is step n+1
  if n==1 
    % Get started with a midpoint method step
    next=pos+tau*vel+0.5*tau^2*accel;
  else
    next=2*pos-prev+tau^2*accel;
  end
  
  % Verlet method: velocity at step n 
  if n>1 
    vel=(next-prev)/(2*tau);
  end
  speed=norm(vel);
  
  % Calculate total energy at step n and store
  energy(n)=0.5*speed^2-1/r(n);
 
  % Update
  prev=pos;
  pos=next;
end

% Plot energy versus time
figure(2);
plot(time,energy);
xlabel('Time (non-dim.)');
ylabel('Total energy (non-dim.)');

% Plot distance (r) versus time
figure (3);
plot(r,time);
xlabel('Time (non-dim.)');
ylabel('Distance (non-dim.)');
