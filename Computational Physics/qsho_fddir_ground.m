% qsho_fddir.m
% Solve the quantum simple harmonic oscillator problem by 
% directly calculating the eigenvalues and eigenvectors of 
% the linear system obtained by finite differencing

% Clear memory and show only a few digits
clear all; format short;

% Define values of the independent variable
h=0.25; 
x=-5:h:5;

% Size of the linear system
L=length(x);

% Construct matrix D2 associated with the derivative
D2=-2*eye(L);
D2=D2+diag(ones(L-1,1),+1)+diag(ones(L-1,1),-1);
D2=D2/h^2;

% Add diagonal component associated with the potential
A=-D2+diag(x.^2);

% Obtain eigenvalues and eigenvectors of A
[V D]=eig(A);

% Analytic Solution
psi = exp(-0.5*x.^2);

% Plot the first four wave functions
figure(1);

  plot(x,V(:,1)/0.3765,'o',x,psi,'x');
  xlabel('x');
  ylabel('\Psi');
  title('Ground state');

% The analytic eigenvalues: odd integers
exact_evals=1:2:2*L-1;

% The numerical estimates
evals=diag(D);