%% Schrodinger's equation for a Gaussian wavepacket with a potential barrier
% Solve the 1-D time-dependent Schrodinger equation for an initial 
% Gaussian wave packet using Crank-Nicolson with a potential barrier

%% Global settings and constants

% Clear memory and show only a few digits
clear all; format short;

% Time step and spatial step
tau=1.e-4;
h=0.0025;

% Parameters of initial wave function
k0=75;  % Average wavenumber
s0=0.05; % Width of Gaussian

% Height of potential barrier
V0=0.25e+4; 

% Total integration time and number of steps
tint=(1+h)/k0;
nsteps=floor(tint/tau)+1;

% Vector of x values
x=0:h:1;
L=length(x);

%% Construct the Hamiltonian

H=-2*eye(L);
H=H+diag(ones(L-1,1),+1)+diag(ones(L-1,1),-1);

% ...with periodic boundary conditions
H(1,L)=1;
H(L,1)=1;

H=-0.5*H/h^2;

% Construct potential matrix corresponding to barrier
V=zeros(L,1);
ii=find((x > 0.6) & (x < 0.7));
V(ii)=V0;
Vmat=diag(V);

% Add potential to Hamiltonian matrix
H=H+Vmat;

% Construct matrix for the linear system solved at each step of 
% Crank-Nicolson. Note that i is sqrt(-1) by default
A=0.5*(eye(L)+0.5*i*tau*H);

%% Define initial conditions and scale

% Initial wave function
C=1./sqrt(s0*sqrt(pi));
psi=C*exp(i*k0*x'); % Oscillatory part
psi=psi.*exp(-0.5*((x-0.25)/s0)'.^2); % Gaussian envelope
psi0=psi;

% Scale for axis
max_psi=max(abs(psi));

%% Crank Nicolson Scheme

time(1)=0;
figure(1);
for n=1:nsteps

    % Update the time
    time(n+1)=time(n)+tau;
    
    % Perform Crank-Nicolson update
    chi=A\psi;
    psi=chi-psi;
    
    % Plot psi versus position and annotate
    subplot(4,1,1),plot(x,real(psi),'b',x,imag(psi),'r');
    axis([0 1 -max_psi max_psi]);
    title(['Blue: Re(\psi),   Red: Im(\psi)      Time: ',...
        num2str(time(n+1))]);
    ylabel('\psi');
    subplot(4,1,2),plot(x,abs(psi).^2,'g')
    axis([0 1 0 max_psi^2]);
    ylabel('|\psi|^2');
    subplot(4,1,3),plot(x,V,'k');
    xlabel('Position');
    ylabel('Potential');
    %3D plot
    subplot(4,1,4),plot3(x,real(psi),imag(psi));
    drawnow;
    
    % Wait for a user response at the start
    if (n == 1) 
      disp('Press a key to continue the calculation...');
      pause;
    end
    
    % Generate a GIF of the plot
    frame = getframe(gcf);
    img =  frame2im(frame);
    [img,cmap] = rgb2ind(img,256);
    if n == 1
        imwrite(img,cmap,'animation.gif','gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(img,cmap,'animation.gif','gif','WriteMode','append','DelayTime',0.1);
    end
    
end