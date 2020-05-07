% qsho_fdover.m
% Solve the quantum simple harmonic oscillator for the lowest
% energy wave function using inverse power iteration relaxation 
% of the linear system obtained by finite differencing, and
% accelerate convergence via over-relaxation

% Clear memory and show only a few digits
clear all; format short;

% Define values of the independent variable
h=0.25;
x=-5:h:5;

% Size of the linear system
L=length(x);

% Number of matrix iterations
nits=10;

% Over-relaxation factor
omega=1;
%omega=1.05;
%omega=1.1;
omega=1.2;
omega=1.5;
%omega=3;

% Construct matrix D2 associated with the derivative
D2=-2*eye(L);
D2=D2+diag(ones(L-1,1),+1)+diag(ones(L-1,1),-1);
D2=D2/h^2;

% Add diagonal component associated with the potential
A=-D2+diag(x.^2);

% Compute inverse matrix
B=inv(A);

% Initial guess for the wave function
%psi0=cos(0.5*pi*x/5)';
%psi0=sin(pi*x/5)';
psi0=zeros(L,1);
psi0(round(2*L/3),1)=1;
psi=psi0;

% Analytic solution for lowest energy wave function
xa=-5:0.1:5;
psi_an=exp(-0.5*xa.^2);

figure(1);
for n=1:nits
  
  % Plot wave function at each iteration, initial wave 
  % function, and analytic solution
  plot(x,psi,'ko',x,psi0,'r*',xa,psi_an)
  xlabel('x'); ylabel('\Psi');
  anno=legend('Numerical solution','Initial guess',...
      'Analytic solution');
  set (anno,'Box','off','Location','NorthWest')
  title(['Iteration: ',num2str(n)]);
  drawnow;
  pause(0.5)
  
  % Inverse iteration (with over-relaxation) & normalisation
  psi_old=psi;
  psi=(1-omega)*psi+omega*B*psi;
  psi=psi/max(abs(psi));
  
  % Determine fractional change
  delta=psi_old-psi;
  fchange(n)=sum(abs(delta./psi));
end

% Determine eigenvalue
energy_est=mean(psi./(B*psi));

disp(['Estimate of eigenvalue E: ',num2str(energy_est)]);

% Plot change versus iteration
figure(2);
semilogy(1:nits,fchange,'ro-');
xlabel('Iteration'); ylabel('Fractional change');
title(['\omega = ',num2str(omega)]);
