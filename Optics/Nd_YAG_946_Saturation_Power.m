clear all
format short

%Estimating the saturation power of a quasi-three-level nd:YAG solid state
%laser operating at 946nm with output coupling R=97%, cavity length l,
%population inversion N, homogenous and uniform gain of g=0.15 (m-1)
%accross the active medium (within the Rayleigh range of the Gaussian beam).

%Globals
R=0.97;
l=0.22;
N=1.8e20;
h=6.626e-34;
c=3e8;
fl=0.0076; %lower lasing fractional population
fu=0.605; %upper lasing fractional population
g=0.15;
tau=230e-6;
Ps=-(1/(8*(fl+fu)*g*tau))*log(R)*N*l*h*c;

C='the saturation power power is:';
string(C);
Ps