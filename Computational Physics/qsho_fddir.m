% qsho_fddir.m
% Solve the quantum simple harmonic oscillator problem by 
% directly calculating the eigenvalues and eigenvectors of 
% the linear system obtained by finite differencing

% Clear memory and show only a few digits
clear all; format short;

% Define values of the independent variable
h=0.0312; 
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

% Plot the first four wave functions
figure(1);
for n=1:4
  subplot(2,2,n),plot(x,V(:,n),'o');
  xlabel('x');
  ylabel('\Psi');
  title(['E = ',num2str(D(n,n))]);
end

% The analytic eigenvalues: odd integers
exact_evals=1:2:2*L-1;

% The numerical estimates
evals=diag(D);

% Plot all the eigenvalues, and compare with the correct values
figure(2);
plot(1:L,evals,'o',1:L,exact_evals,'r+');
xlabel('Energy state number'); ylabel('Energy eigenvalue');
anno=legend('Finite difference result','Exact eigenvalue');
set (anno,'Box','off','Location','SouthEast')