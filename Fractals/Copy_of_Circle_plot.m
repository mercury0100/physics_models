clear all, format short

%define angular resolution and sigma and rho values
res=30;
sigvals=0:(1/res):pi;
rhovals=0:(1/res):pi;

%define radius
r=1;

figure(1);

for m=1:res

for n=1:res
 
%define spherical coordinates and compression factor
k=rhovals(m)/sigvals(n);

x(n,m)=k*r*sin(sigvals(n)*2*pi)*cos(rhovals(m)*2*pi);
y(n,m)=k*r*sin(sigvals(n)*2*pi)*sin(rhovals(m)*2*pi);
z(n,m)=r*cos(sigvals(n)*2*pi);

plot3(x,y,z);
drawnow;

end

end
