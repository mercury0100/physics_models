% Clear memory and show only a few digits
clear all; format short;

% advection_lw.m
% Solve the advection equation using Lax-Wendroff with periodic BCs

% Thermal conductivity and c
c=1;
kappa=0.1;

% Spatial and time steps
tau=0.002;


% Time step, total time and spatial step
tau=2e-3;
L=20;
h=1/(L-1);
tf=0.5*(1+h)/c;

% Initial profile constant sigma
sig = 0.1;

% Number of time steps
nsteps=floor(tf/tau);

% Vector of x values
jvals=1:L;
x=(jvals-1)/h;

% Diagonal matrix to implement the centred spatial first-derivative
D=diag(ones(L-1,1),+1)-diag(ones(L-1,1),-1);

% define g and f
g=0:0.01:1;
f=0:0.01:1;
gnsteps=length(g);
fnsteps=length(f);

% Initial conditions (a Gaussian pulse at x = 0.5)
amp=exp(-0.5*(x-0.5).^2/sig^2)';
amp0=amp; 

% Record a(x,t) matrix for visualisation
amp_xt(:,1)=amp0;

% Additional elements for periodic boundary conditions
D(1,L)=-1;
D(L,1)=+1;

% Lax-Wendroff second derivative matrix
D2=abs(D)-2*eye(L);

for fn=1:fnsteps
    
for gn=1:gnsteps

% Update matrix
A=eye(L)-0.5*g(gn)*D+0.5*f(fn)*D2;

% Calculate the spectral radius of M
rho(fn,gn)=max(abs(eig(A)));

end
    
end

% Visualisation of amplitude versus position and time. MATLAB
% uses a surfc(x(i),y(j),z(j,i)) convention, hence the transpose
figure(1);
colormap(summer);
surfc(g,f,rho');
shading interp;
xlabel('g'); ylabel('f'); zlabel('Spectral Radius');
