clear all, format short

%define angular resolution and sigma values
res=500;
sigvals=0:(1/res):2*pi;

%define radius
r=1;

figure(1);
for n=1:res
    
%define circle coordinates
x(n)=r*cos(sigvals(n)*2*pi);
y(n)=r*sin(sigvals(n)*2*pi);

plot(x,y,'o');
drawnow;

end

