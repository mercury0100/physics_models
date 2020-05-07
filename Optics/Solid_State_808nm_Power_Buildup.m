clear all
format short

%Modeling the power buildup of a nd:YAG solid state laser operating at
%946nm with uniform active gain of g (m-1), logarithmic cavity loss gamma
%(m-1), optical coupler efficiencies R1 and R2 and initial power Pi

%global variables
go=1e4;
lg=280e-6;
la=350e-6;
gamma=3.5e3;
R1=0.99;
R2=0.33;
tau=230e-6;
h=6.626e-34;
v=3e8/(808e-9);
P(1)=h*v/tau;
Ps=0.010;

%Initiate loop
n=1;

%Calculate Cavity power as a function of roundtrips

while n<200
    g(n+1)=go/(1+2*(P(n)/Ps));
    P(n+1)=P(n)*(R2*R1*exp((g(n+1)*lg-gamma*la)*2));
    Pout(n+1)=P(n)*(1-R2)*(R1*exp((g(n+1)*lg-gamma*la)*2));
    k(n+1)=n;
    n=n+1;
end

figure(1)
subplot(3,1,1)
plot(k,P,'r')
title('Cavity Power (W)')
subplot(3,1,2)
plot(k,Pout,'b')
title('Output Power (W)')
subplot(3,1,3)
plot(k,g,'g')
title('Gain')