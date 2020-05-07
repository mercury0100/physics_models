% poisson_sor.m
% Apply successive overrelaxation to the Dirichlet problem for the 
% Poisson equation

% Clear memory and show only a few digits
clear all; format short; clf;

% Parameters
h=0.025;         % Spatial step
max_iter=1e+4;   % Maximum number of iterations
min_diff=1e-6;   % Convergence criterion
skip=1;         % Only plot every skip iterations

% Vectors of x and y values
x=0:h:1; y=x;
L=length(x);

% Analytic result for optimal value
omega_opt=2/(1+sin(pi/L));
disp('Optimal omega:');
disp(omega_opt)
omega=omega_opt;

% Initial phi
phi=zeros(L);

% Charge density: a delta function at the centre of the region
sigma=zeros(L);
sigma(round(L/2),round(L/2))=1/h^2;

% Relaxation loop
figure(1);
for iter=1:max_iter
    
  % Keep this for the convergence test
  old_phi=phi;
  
  % Update in place using SOR, preserve BCs
  for j=2:L-1
    for l=2:L-1
      gs=0.25*(phi(j-1,l)+phi(j+1,l)+ ...
                     phi(j,l-1)+phi(j,l+1)+h^2*sigma(j,l));
      phi(j,l)=(1-omega)*phi(j,l)+omega*gs;
    end
  end

  % Plot updated solution
  if rem(iter,skip) == 0
    surfc(x,y,phi');
    %shading interp;
    xlabel('x');
    ylabel('y');
    zlabel('Potential \phi');
    title(['Iteration: ',num2str(iter)]);   
    drawnow;
  end
  
  % Break if change is suitably small
  diff=max(abs(phi-old_phi));
  if (diff < min_diff) 
    disp('Relaxation converged...');
    break;
  end
end

% Display number of iterations
disp('Number of iterations: ');
disp(iter);

% Determine the electric field by differencing. Note MATLAB's
% perverse identification of the first index in a matrix with
% y, wherever a functional dependence z=z(x,y) is implied.
[Ey, Ex]=gradient(phi,h);
Ex=-Ex; Ey=-Ey;

% Calculate the (approximate) analytic solution
Exa=zeros(L);
Eya=zeros(L);
for j=1:L
  for l=1:L
    r=sqrt((x(j)-0.5)^2+(y(l)-0.5)^2);
    fac=1/(2*pi);
    if (r ~= 0)
      Exa(j,l)=fac*(x(j)-0.5)/r^2;
      Eya(j,l)=fac*(y(l)-0.5)/r^2;
    end
  end
end

% Compare the electric field obtained by the SOR solution with
% the (approximate) analytic solution
figure(2); clf reset;
hold on;
contourf(x,y,phi');
scale=2;
quiver(x,y,Ex',Ey',scale,'w');
quiver(x,y,Exa',Eya',scale,'r');
xlabel('x');
ylabel('y');
title('Electric field {\bf E} and contours of \phi');
hold off;