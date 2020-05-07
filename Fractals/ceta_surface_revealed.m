% n iterations of x^2+c
clear all; format rat; format short

%define c, boundaries and number of iterations
c_max=1000;
cstepsize=1;   
cvals=-c_max:cstepsize:c_max;
csteps=length(cvals);
nsteps=7;
nvals=0:1:nsteps;

for m=1:csteps-1

%define initial entry x and "it" vector
x=0;
it(1,1)=x;

for n=1:nsteps
    
    x=x^2+cvals(m);
    
    it(m+1,n+1)=x;
    
end

end



% Visualisation of iterations
figure(1);
colormap(summer);
surfc(cvals,nvals,it');
%set(gca, 'zscale', 'log');
shading interp;
xlabel('c'); ylabel('n'); zlabel('cprimes');


    