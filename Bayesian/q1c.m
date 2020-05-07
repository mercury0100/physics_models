clear all
close all
format short

vo = 1.65;
g = 9.8*10^-3;
Data = [51,198,147,39,144,288,226,13,89,137,98,106,102,100,55,44,129,22,77,82];
[X,Y] = meshgrid(50:0.1:250,0.05:0.001:0.5);
L  = ones;

for i=1:length(Data)
    r = Data(i);
    v = sqrt(g*(r+X));
    Upper = g*exp((-(v - vo).^2)./(2*Y.^2));
    Lower = 2*Y.*sqrt(2*pi*g*(r+X));
    likelihood = Upper./Lower;
    L = L.*likelihood;
end 

contour(X,Y,L,20)

