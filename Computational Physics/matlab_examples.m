% matlab_examples.m
% 
% An example of a script: the MATLAB examples from lecture 1. This 
% script is not intended to be run independently (although it will
% run without error) it is intended to collect the examples.

% Clear memory and only show a few digits
clear all; format short;

% Data

x=1;
whos x;

y=[1 2 3]; z=y';
disp(y); disp(z);

% Array operations, constants, digits of precision
y^2 % A common error
y.^2
z*y
v=1:10

w=pi
format long 
w
format short

% Use of ellipsis to line break
x=...
    2;

% Graphical help browser
doc

% Program control

% If-else statements
 x=input('Value of integer x? ');
if x==0 
  disp('x is equal to zero');
else
  if x>0 
    disp('x is positive');
  else 
    disp('x is negative');
  end
end

% Loops
tau=0.1; LENGTH=50000;
for n=1:LENGTH 
  time(n)=tau*n;
end

% Array version - use this!
time2=tau*(1:LENGTH);

% Find elements of an array matching a condition
find((time > 1) & (time < 1.2))

% Example in lecture
x1=-0.5:0.05:0.5; x2=x1;
L=length(x1);
x2(find(abs(x1)<0.25))=0;
plot(1:L,x1,'o',1:L,x2,'+');
xlabel('Array index');
ylabel('Array value');

% Matrix algebra
A=rand(3,3) % This is the same as rand(3)
b=rand(3,1)
x=A\b
res=A*x-b
max(res)

% Compare with the machine precision/epsilon eps...
disp(['eps = ',num2str(eps)]);

% What is the machine precision?
disp(['1+0.5*eps-1 = ',num2str(1+0.5*eps-1)]);

% Example of calculating determinant, eigenvectors and values
det(A)
eig(A)
[V,D]=eig(A)

% 1-D plotting example
x=0:0.1:2*pi; f1=sin(x); f2=cos(x);
plot(x,f1,x,f2);
xlabel('x')
ylabel('y')
title('A simple plot in MATLAB... always label your axes!');

% 2-D plotting examples
colormap(hot); % Avoid the rainbow colormap!
y=-pi:0.1:pi;
z=sin(x)'*cos(y);

% Filled contours...
contourf(x,y,z);
axis([min(x) max(x) min(y)...
   max(y)]);
xlabel('x'); ylabel('y');

% Or a surface
surfc(x,y,z);
axis([min(x) max(x) min(y) max(y)]);
xlabel('x'); ylabel('y'); zlabel('z');

% Functions - see the file dice.m
N=10000;
test=dice(N);     % Generate 10000 throws of a die
x=1:6;            % Horizontal axis values
hist(test,x)      % Histogram of results
xlabel('Value'); 
ylabel('Number');

% Defining and calling an anonymous function
x=[0:0.1:10];
cauchy=@(x,a) 1./pi./(1+(x-a).^2);
plot(x,cauchy(x,5));
xlabel('x');
ylabel('y');








