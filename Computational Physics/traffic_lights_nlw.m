% traffic_lights_nlw.m
% Solve the fluid model describing traffic for a line of cars 
% starting off at traffic lights, using nonlinear Lax-Wendroff 
% with periodic BCs

% Clear memory and show only a few digits
clear all; format short;

% Spatial step
h=0.01;

% Generalized CFL condition (tau/h=1, since max(c)=1)
tau=h; 
fac=tau/h;

% Perform two loops of the periodic domain, assuming a maximum speed
nloops=2;
nsteps=nloops*ceil((1+h)/tau); 

% Vector of x values
x=0:h:1;
L=length(x);

% Diagonal matrix to implement the spatial first-derivative
D=diag(ones(L-1,1),+1)-diag(ones(L-1,1),-1);

% Additional elements for periodic boundary conditions
D(1,L)=-1;
D(L,1)=+1;

% Lax-Wendroff second derivative matrix
D2=abs(D)-2*eye(L);

% Initial conditions on density
eps=0.05;
rho=zeros(L,1);
ii=find((x > 0.25) & (x <= 0.5-eps));
rho(ii)=1;
ii=find((x > 0.5-eps) & (x <= 0.5));
rho(ii)=linspace(1,0,length(ii));

% Initial wave speed and flux profiles
c=1-2*rho;
F=rho.*(1-rho);

% Record wave amplitude for visualisation
rho_xt(:,1)=rho;

% March forwards in time
figure(1);
time(1)=0;
for n=1:nsteps
    
    % Nonlinear Lax-Wendroff method
    term1=-0.5*fac*D*F;
    term2=1-0.25*fac*D*c;
    term3=0.5*fac^2*D2*F;

    rho=rho+term1.*term2+c.*term3;
    c=1-2*rho;
    F=rho.*(1-rho);
    
    % Update the time
    time(n+1)=time(n)+tau;

    % Plot the density profile
    plot(x,rho,'ro-');
    xlabel('Position');
    ylabel('Density');
    title(['Time: ',num2str(time(n+1))]);      
    axis([0 1 0 1]);
    pause(0.1);
    drawnow;
    
    % Record density for visualisation
    rho_xt(:,n+1)=rho;
end

% Visualisation of amplitude versus position and time. MATLAB
% uses a contourf(x(i),y(j),z(j,i)) convention, hence the transpose
% Visualisation of density versus position and time
figure(2);
colormap(summer);
contourf(x,time,rho_xt');
xlabel('Position');
ylabel('Time');
zlabel('density');
xlabel('Position'); ylabel('Time'); zlabel('Amplitude');
