% n iterations of x^2+c
clear all; format rat;

%define c and number of iterations
c=3;   
nsteps=50;




%define initial entry x and vector containing successive iterations
x=0;
primes(1)=x;

for n=1:nsteps;
    
    x=x^2+c;
    
    primes(n+1)=x;
    
end

primes



    