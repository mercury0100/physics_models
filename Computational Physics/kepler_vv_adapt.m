% kepler_vv_adapt.m
% Motion under a central force, using the modified Verlet method
% (Velocity-Verlet), with the Numerical Recipes step doubling 
% method.

% Clear memory and only show a few digits
clear all; format short;

% Nominal time step (non-dim.)
tau00=0.1;

% Adaptive timestep parameters
DEL0=1.e-10;
min_tau=1.e-16;
max_tau=1000;
k=3; % Order of Velocity-Verlet in one step

% Display parameter
fspec = '%7.3e';

% Initial position (non-dim.) - this should be fixed
pos=[1 0];

% Initial velocity (non-dim.) - vary the y-component
vel=[0 1.4];

% Total integration time
T=800;

% Number of integration steps
steps=ceil(T/tau00);

% Avoid plotting every step
skip=ceil(steps/10);
                                 
% Calculate trajectory from analytic solution. See Appendix
% of Week 2 lecture 1.
[xan,yan]=kepler_analytic(vel,T);

% Modified Verlet integration
figure(1);
time(1)=0;
energy(1)=0.5*norm(vel)^2-1/norm(pos);
n=1; % Step counter
while time < T
  
  % Always start with default time step
  tau=tau00;

  x(n)=pos(1);
  y(n)=pos(2);

  % Plot positions every so many steps
  if (rem(n,skip)==0)
    plot(x,y,'g-',pos(1),pos(2),'ko',xan,yan,'b',0,0,'ro')
    title(['Time: ',num2str(time(n),fspec)])
    xlabel('x');
    ylabel('y');
    axis equal; % Preserve aspect ratio
    drawnow;
  end

  % Step doubling loop: take a full step then two half steps
  ACC=0;
  m=0; % Counter of number of times step estimate is made
  while (ACC == 0)

    % Velocity-Verlet method as one step
    accel=-pos/norm(pos)^3;
    pos_next1=pos+tau*vel+0.5*tau^2*accel;
    accel_next=-pos_next1/norm(pos_next1)^3;
    vel_next1=vel+0.5*tau*(accel+accel_next);

    % Velocity-Verlet method as two half steps
    tau2=0.5*tau;
    pos_next21=pos+tau2*vel+0.5*tau2^2*accel;
    accel_next21=-pos_next21/norm(pos_next21)^3;
    vel_next21=vel+0.5*tau2*(accel+accel_next21);

    pos_next22=pos_next21+tau2*vel_next21+0.5*tau2^2*accel_next21;
    accel_next22=-pos_next22/norm(pos_next22)^3;
    vel_next22=vel_next21+0.5*tau2*(accel_next21+accel_next22);

    % Difference in final positions
    delta=norm(pos_next1-pos_next22);

    % Exit step doubling loop, or calculate estimate for an accurate 
    % step
    if (delta <= 1.5*DEL0)
      % Accept and make updates
      ACC=1;
      pos=pos_next22;
      vel=vel_next22;
      taulog(n)=tau;
      n=n+1;
      time(n)=time(n-1)+tau;
      energy(n)=0.5*norm(vel)^2-1/norm(pos);
    else 
      % Calculate step size to achieve accuracy criterion
      tau0=tau*(DEL0/delta)^(1/k);
      tau=0.5^m*tau0;
      % disp(['Reducing time step to = ',num2str(tau)]);
      if (m > 0) 
        disp(['Using m > 0. Criterion too stict?']); 
      end
      m=m+1;
    end

    % Terminate if tau too small/large
    if ((tau < min_tau) | (tau > max_tau))
      disp('Timestep exceeded allowed range...');
      break
    end

  end
end

% Plot energy and of the number of times the step is halved vs. time
figure(2);
subplot(2,1,1);
plot(time,energy);
axis([0 T min(energy) max(energy)]);
ylabel('Total energy (non-dim.)');
subplot(2,1,2);
semilogy(time(1:end-1),taulog);
xlabel('Time (non-dim.)');
ylabel('Time step (non-dim.)');
axis([0 T min(taulog) max(taulog)]);
