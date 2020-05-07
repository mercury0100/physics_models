% diffusion_ftcs.m
% Solve the 1-D diffusion equation for an initial spike profile 
% with Dirichlet conditions using FTCS, in a matrix formulation

% Clear memory and show only a few digits
clear all; format short;

% Thermal conductivity
kappa=1;

% Time step and spatial step
tau=1e-4;
h=0.02;

% Number of time steps
nsteps=25;

% Calculate the ratio tau/(th/2), where th is the approximate diffusion
% time for one spatial step h, and display
th=h^2/kappa;
disp(['Ratio tau/(0.5*th): ',num2str(tau/(0.5*th))]);

% Vector of x values
x=0:h:1;
L=length(x);

% Construct the matrix D associated with the second spatial 
% derivative and the boundary conditions
D=-2*eye(L);
D=D+diag(ones(L-1,1),+1)+diag(ones(L-1,1),-1);
D=kappa*tau*D/h^2;

% Impose the Dirichlet boundary conditions
D(1,:)=zeros(1,L);
D(L,:)=zeros(1,L);

% Construct the update matrix
A=eye(L)+D;

% Initial conditions: a spike at x = 1/2
temp1=zeros(L,1);
temp1(round(L/2))=1/h;
temp=temp1;

% Record T(x,t) matrix for visualisation
temp_xt(:,1)=temp;

% March forwards in time
time(1)=0;
figure(1);
for n=1:nsteps

    % Update the time and the temperature profile
    time(n+1)=time(n)+tau;
    temp=A*temp;
    
    % Calculate the profile for the (approximate) analytic solution
    sig=sqrt(2*kappa*time(n+1));
    temp_an=exp(-(x-0.5).^2/(2*sig^2))/(sqrt(2*pi)*sig);
    
    % Plot the temperature versus position, the initial temperature
    % profile, and the analytic values, with annotations
    plot(x,temp,'ro-',x,temp1,'g-*',x,temp_an)
    xlabel('Position (non-dim.)');
    ylabel('Temperature (non-dim.)');
    handle=legend('Numerical solution','Initial profile',...
      'Analytic solution');
    set (handle,'Box','off','Location','NorthWest')
    title(['Time: ',num2str( time(n))]);
    if (n > nsteps/2) % Follow evolution more closely near end
      axis([0 1 min(temp_an) max(temp_an)]);
    end
    drawnow; pause(0.01);
    
    % Record T(x,t) matrix for visualisation
    temp_xt(:,n+1)=temp;
    
end

% Visualisation of temperature versus position and time. MATLAB
% uses a surfc(x(i),y(j),z(j,i)) convention, hence the transpose
figure(2);
surfc(x,time,temp_xt');
%shading interp;
xlabel('Position'); ylabel('Time'); zlabel('Temperature');
