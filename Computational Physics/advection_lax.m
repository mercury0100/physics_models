% advection_lax.m
% Solve the advection equation using the Lax method with periodic BCs

% Clear memory and show only a few digits
clear all; format short;

% Propagation speed
c=1;

% Spatial and time steps
h=0.01;
tau=0.005;

% Number of steps such that the pulse should propagate through
% the periodic domain once
nsteps=ceil((1+h)/(c*tau));

% Vector of x values
x=0:h:1;
L=length(x);

% Diagonal matrix to implement the centred spatial first-derivative
D=diag(ones(L-1,1),+1)-diag(ones(L-1,1),-1);

% Additional elements for periodic boundary conditions
D(1,L)=-1;
D(L,1)=+1;

% Lax averaging matrix
A=0.5*abs(D);

% Update matrix
M=A-0.5*(c*tau/h)*D;

% Initial conditions (a Gaussian pulse at x = 0.5)
sig=0.1;
amp=exp(-0.5*(x-0.5).^2/sig^2)';
amp0=amp; 

% Record a(x,t) matrix for visualisation
amp_xt(:,1)=amp0;

% March forwards in time
figure(1);
time(1)=0;
for n=1:nsteps
    
    % Update the time and the wave amplitude profile
    time(n+1)=time(n)+tau;
    amp=M*amp;

    % Calculate the profile for the exact analytic solution
    k=floor(c*time(n+1)/h); % The index the peak has reached
    xp=x([k+1:end 1:k]); % Shift x to left by k steps
    
    % Plot the amplitude versus position, the initial amplitude
    % profile, and the analytic values, with annotations
    plot(x,amp,'bo',x,amp0,'r',xp,amp0,'g+')
    xlabel('Position (non-dim.)');
    ylabel('Amplitude (non-dim.)');
    title(['\tau/\tau_{\rm max}= ',num2str(tau/(h/c)),...
        '     Time: ',num2str(time(n))]);
    drawnow; 
    pause(0.01);
    
    % Record a(x,t) matrix for visualisation
    amp_xt(:,n+1)=amp;
    
end

% Add annotation
handle=legend('Numerical solution','Initial profile',...
  'Analytic solution');
set (handle,'Box','off','Location','NorthWest')

% Visualisation of amplitude versus position and time. MATLAB
% uses a surfc(x(i),y(j),z(j,i)) convention, hence the transpose
figure(2);
colormap(summer);
surfc(x,time,amp_xt');
shading interp;
xlabel('Position'); ylabel('Time'); zlabel('Amplitude');

% Calculate the spectral radius of M
rho=max(abs(eig(M)));
disp(['Spectral radius: ',num2str(rho)]);
