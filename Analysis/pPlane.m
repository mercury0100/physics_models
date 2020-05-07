clear all; format short; hold on;
[X,Z] = meshgrid(-5:.2:5,-5:.2:5);
dX = X*(1-X)-0.2*Z;
dZ = (X-Z)/25;
quiver(X,Z,dX,dZ)
X = -1:.05:1; nullX = (X-X.^2)/0.2; nullZ = X; plot(X,nullX,X,nullZ)
% phase diagram with vectorfield