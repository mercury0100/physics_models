% n iterations of x^2+c
clear all; format rat; format short

%define c and number of iterations
csteps=5;   
nsteps=5;
ceta=0:1:nsteps;
neta=0:1:nsteps;

for c=1:csteps

%define initial entry x and vector containing successive iterations
x=0;
it(1,1)=x;

for n=1:nsteps
    
    x=x^2+c;
    
    it(c+1,n+1)=x;
    
end

end

% Visualisation of iterations
figure(1);
colormap(summer);
surfc(ceta,neta,it);
set(gca, 'zscale', 'log');
shading interp;
xlabel('ceta'); ylabel('n'); zlabel('cprimes');


    