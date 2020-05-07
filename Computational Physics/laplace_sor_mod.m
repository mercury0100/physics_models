% laplace_sor.m
% Solve the Laplace equation with Dirichlet BCs using successive
% overrelaxation (SOR)

% Clear memory and show only a few digits
clear all; format short; clf; 

% Parameters
h=0.025;         % Spatial step
max_iter=1e+4;  % Maximum number of iterations
min_diff=1e-6;  % Convergence criterion
skip=10;        % Only plot every skip iterations

hvals=[0.025,0.05,0.1,0.2];
itern=length(hvals);
for n=1:itern;
    h=hvals(n);

% Vectors of x and y values
x=0:h:1; y=x;
L=length(x);

% Analytic result for optimal value
omega_opt=2/(1+sin(pi/L));
disp(['Optimal omega:',num2str(omega_opt)]);
omega=omega_opt;

% Over-relaxation factor
%omega=1.74;
%omega=1.2;

% Initial phi
phi=zeros(L);

% Impose BCs
phi(2:L-1,L)=1;

% Relaxation loop
figure(1);
colormap(summer);
for iter=1:max_iter
    
  % Keep this for the convergence test
  old_phi=phi;
  
  % Loop over interior points: update in place using SOR, preserve BCs
  for j=2:L-1
    for l=2:L-1
      gs=0.25*(phi(j-1,l)+phi(j+1,l)+ ...
                     phi(j,l-1)+phi(j,l+1));
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

figure(3)
Lvals(n)=L;
itervals(n)=iter;
loglog(Lvals,itervals,'o');

end

% Display number of iterations
disp(['Number of iterations: ',num2str(iter)]);

% Determine the electric field by differencing. Note MATLAB's
% perverse identification of the first index in a matrix with
% y, wherever a functional dependence z=z(x,y) is implied.
[Ey, Ex]=gradient(phi,h);
Ex=-Ex; Ey=-Ey;

% Plot also the electric field 
figure(2); clf reset;
hold on;
colormap(summer);
contourf(x,y,phi');
scale=10;
quiver(x,y,Ex',Ey',scale,'w');
xlabel('x');
ylabel('y');
title('Electric field {\bf E} and contours of \phi');
hold off;
