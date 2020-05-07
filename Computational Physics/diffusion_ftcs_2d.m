% diffusion_ftcs_2d.m
% Solve the 2-D diffusion equation for an initial heat spike with 
% Dirichlet conditions using FTCS

% Clear memory and show only a few digits
clear all; format short;

% Thermal conductivity
kappa=1;

% Time step and spatial step
tau=0.5e-4;
h=0.02;

% Number of time steps
nsteps=100; 

% Check stability criterion
fac=4*kappa*tau/h^2;
disp(['Stability factor: ',num2str(fac)]);

% Values of x and y
x=0:h:1; y=0:h:1;
L=length(x);

% Initial conditions: a spike at x = y = 1/2
temp=zeros(L,L);
temp(round(L/2),round(L/2))=1/h^2;

% March forwards in time
time=0;
figure(1);
for n=1:nsteps
    
    % Update the time
    time=time+tau;
    
    % Update the T(x,y) values with a nested loop
    next=temp;
    for l=2:L-1
      for m=2:L-1
        terms=temp(l-1,m)+temp(l+1,m)+temp(l,m-1)+...
            temp(l,m+1)-4*temp(l,m);
        next(l,m)=temp(l,m)+kappa*tau*terms/h^2;
      end
    end
    temp=next;
    
    % Surface plot of T(x,y) with annotations
    surfc(x,y,temp');
    %shading interp;
    title(['Time: ',num2str(time)]);
    xlabel('x'); ylabel('y'); zlabel('Temperature');
    if (n < nsteps/2)
      axis([0 1 0 1 0 0.1*1/h^2]);
    else % Follow evolution more closely near end
      axis([0 1 0 1 0 max(max((temp)))]);
    end
    drawnow;

end
