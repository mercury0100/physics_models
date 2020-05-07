% poisson_fft.m
% Solve the Poisson equation using FFTs (and hence implicit periodic
% boundary conditions)

% Clear memory and show only a few digits
clear all; format short;

% Spatial step
h=0.025;

% Vectors of x and y values
x=0:h:1; y=x;
L=length(x);

% Initial phi
phi=zeros(L);

% Charge density: a delta function at the centre of the region
sigma=zeros(L);
sigma(round(L/2),round(L/2))=1/h^2;
%sigma(L,L)=1/h^2;

% Fourier transform charge density
fsigma=fft2(sigma);

% Construct arrays for FTs of phi, Ex, Ey
fphi=zeros(L);
fEx=zeros(L);
fEy=zeros(L);
  
% Loop over m and n and construct FTs of phi, Ex, Ey
for m=1:L
  for n=1:L
      
    % Denominator in FT of solution at grid point
    den=cos(2*pi*(m-1)/L)+cos(2*pi*(n-1)/L)-2;
    if (den ~= 0)
      
      % Construct FT of phi at grid point
      fphi(m,n)=-0.5*h^2*fsigma(m,n)/den;
      
      % Construct FTs of electric field components
      fEx(m,n)=-sqrt(-1)*sin(2*pi*(m-1)/L)*fphi(m,n)/h;
      fEy(m,n)=-sqrt(-1)*sin(2*pi*(n-1)/L)*fphi(m,n)/h;
    
    end
  end
end

% Invert FTs to obtain phi(x,y), Ex(x,y), Ey(x,y)
phi=real(ifft2(fphi));
Ex=real(ifft2(fEx));
Ey=real(ifft2(fEy));

% Show phi
figure(1);
colormap(summer);
surfc(x,y,phi');
%shading interp;
xlabel('x');
ylabel('y');
zlabel('Potential \phi');

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

% Compare the electric field obtained by the Fourier solution with
% the (approximate) analytic solution
figure(2); clf reset;
hold on;
colormap(summer);
contourf(x,y,phi');
scale=2;
quiver(x,y,Ex',Ey',scale,'w');
quiver(x,y,Exa',Eya',scale,'r');
xlabel('x');
ylabel('y');
title('Electric field {\bf E} and contours of \phi');
hold off;
