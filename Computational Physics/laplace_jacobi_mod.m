% laplace_jacobi.m
% Solve the Laplace equation with Dirichlet BCs using Jacobi 
% relaxation

% Clear memory and show only a few digits
clear all; format short; clf;

% Parameters
h=0.1;         % Spatial step
max_iter=1e+4;  % Maximum number of iterations
min_diff=1e-6;  % Convergence criterion
skip=10;        % Only plot every skip iterations


hvals=[0.05];
itern=length(hvals);

for n=1:itern;
    h=hvals(n);
    
    x=0:h:1; y=x;
    L=length(x);
    
    % Initial phi
%phi_new=0.2*ones(L);
phi_new=rand(L);

% Impose BCs
phi_new(2:L-1,L)=1;
phi_new(:,1)=0;
phi_new(1,:)=0;
phi_new(L,:)=0;

%phi_new=rand(L);

% Relaxation loop
figure(1);
colormap(summer);
for iter=1:max_iter
    
  phi=phi_new;
  
  % Loop over interior points and update. Note the boundary
  % values are preserved.
  for j=2:L-1
    for l=2:L-1
      phi_new(j,l)=0.25*(phi(j-1,l)+phi(j+1,l)+ ...
                           phi(j,l-1)+phi(j,l+1));
    end
  end

  % Plot updated solution
  if rem(iter,skip) == 0
    surfc(x,y,phi_new');
    %shading interp;
    xlabel('x');
    ylabel('y');
    zlabel('Potential \phi');
    title(['Iteration: ',num2str(iter)]);   
    drawnow;
  end
  
  % Break if change is suitably small
  diff=max(abs(phi_new-phi));
  if (diff < min_diff) 
    disp('Relaxation converged...');
    break;
  end
end


    
    figure(3);
    Lvals(n)=L;
    itervals(n)=iter;
    loglog(itervals,Lvals,'o');
    
    n=n+1;

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
contourf(x,y,phi');
scale=10;
colormap(summer);
quiver(x,y,Ex',Ey',scale,'w');
xlabel('x');
ylabel('y');
title('Electric field {\bf E} and contours of \phi');
hold off;
