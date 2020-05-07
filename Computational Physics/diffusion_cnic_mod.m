% diffusion_cnic.m
% Solve the diffusion equation using Crank-Nicolson

% Clear memory and show only a few digits
clear all; format short;

% Parameters
kappa=1;          % Thermal conductivity
tau=1e-3;         % Time step
hvals=[0.005,0.01,0.025,0.05,0.1];           % Spatial step
tint=0.5;  % Total integration time
nsteps=floor(tint/tau);       % Number of time steps
skip=1;           % Only plot every now and then

 for m=1:length(hvals);
        
h=hvals(m);
fac=kappa*tau/h^2;

% Display value of FTCS stability factor
disp(['FTCS stability factor: ',num2str(fac)]);

% Vector of x values
x=0:h:1;
L=length(x);

% Construct the matrix D associated with the second spatial derivative
D=-2*eye(L);
D=D+diag(ones(L-1,1),+1)+diag(ones(L-1,1),-1);

% Impose Dirichlet boundary conditions
D(1,:)=zeros(1,L);
D(L,:)=zeros(1,L);

% Construct the matrix for the linear system solved at each step of 
% Crank-Nicolson
A=0.5*(eye(L)-0.5*fac*D);

% Initial conditions
temp=sin(pi.*x)';

% Record T(x,t) matrix for visualisation
%temp_xt(:,1)=temp;

% Analytic Solution


% March forwards in time
time(1)=0;
%figure(1);


 for n=1:nsteps

    % Update the time
    time(n+1)=time(n)+tau;
    
    % Perform Crank-Nicolson update
    chi=A\temp;
    temp=chi-temp;
    
    % Calculate the profile for the (approximate) analytic solution
    sig=sqrt(2*kappa*time(n+1));
    temp_an=exp(-(pi^2)*0.5)*sin(pi.*x);
    
    % Plot the temperature versus position and the analytic values, 
    % with annotations
    %if rem(n,skip) == 0
    %plot(x,temp,'ro-',x,temp_an)
    %xlabel('Position (non-dim.)');
    %ylabel('Temperature (non-dim.)');
    %handle=legend('Numerical solution','Analytic solution');
    %set (handle,'Box','off','Location','NorthWest')
    %title(['f = \kappa\tau/h^2 = ',num2str(kappa*tau/h^2),'     '...
    %    'Time: ',num2str(time(n))]);
      %pause(0.1);
    %drawnow;
    %end
   
    
    % Record T(x,t) matrix for visualisation
    %temp_xt(:,n+1)=temp;
    
    end

    % Error in Final Temp
Er(m)=max(abs(temp-(temp_an')))
    
    end

% Visualisation of temperature versus position and time. MATLAB
% uses a surfc(x(i),y(j),z(j,i)) convention, hence the transpose
%figure(2);
%colormap(summer);
%surfc(x,time,temp_xt');
%shading interp;
%xlabel('Position'); ylabel('Time'); zlabel('Temperature');

figure(3);
loglog(hvals,Er);
