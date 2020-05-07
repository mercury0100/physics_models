% diffusion_fimp.m
% Solve the 1-D diffusion equation using the fully implicit scheme
% with Dirichlet BCs

% Clear memory and show only a few digits
clear all; format short;

% Parameters
kappa=1;           % Thermal conductivity
h=0.02;            % Spatial step
fac=20;            % For FTCS, this must be <= 1/2 for stability
tau=fac*h^2/kappa; % Time step
nsteps=100;        % Number of time steps
tint=tau*nsteps;   % Total integration time
skip=5;            % Only plot every now and then

% Display value of FTCS stability factor
disp(['FTCS stability factor: ',num2str(fac)]);

% Vector of x values
x=0:h:1;
L=length(x);

% Construct the matrix D associated with the second spatial derivative
D=-2*eye(L);
D=D+diag(ones(L-1,1),+1)+diag(ones(L-1,1),-1);

% Construct the matrix in the linear system being solved at each step
M=eye(L)-fac*D;

% Impose Dirichlet boundary conditions
M(1,:)=zeros(1,L);
M(1,1)=1;
M(L,:)=zeros(1,L);
M(L,L)=1;

% Initial conditions and BCs
temp=zeros(L,1);
temp(1)=0;
temp(L)=1;

% Record T(x,t) matrix for visualisation
temp_xt(:,1)=temp;

% March forwards in time
time(1)=0;
figure(1);
%hold on;
for n=1:nsteps

    % Update the time
    time(n+1)=time(n)+tau;   
    
    % Solve the matrix equation for a time step
    temp=M\temp;     
    
    % Plot the temperature versus position and the equilibrium 
    % solution, with annotations
    if (rem(n,skip) == 0) 
      plot(x,temp,'ro',x,x)
      axis([0 1 0 1]);
      xlabel('Position (non-dim.)');
      ylabel('Temperature (non-dim.)');
      handle=legend('Numerical solution','Equilibrium');
      set (handle,'Box','off','Location','NorthWest')
      title(['Time: ',num2str(time(n))]);
      drawnow;
    end
    
    % Record T(x,t) matrix for visualisation
    temp_xt(:,n+1)=temp;
    
end

% Visualisation of temperature versus position and time. MATLAB
% uses a surfc(x(i),y(j),z(j,i)) convention, hence the transpose
figure(2);
colormap(summer);
surfc(x,time,temp_xt');
xlabel('Position'); ylabel('Time'); zlabel('Temperature');
