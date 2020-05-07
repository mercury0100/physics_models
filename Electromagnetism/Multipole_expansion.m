clear all;
format short;

% Calculate and compare the exact electrostatic potential of a water
% molecule against its multipole expansion, with dipole aligned along the z
% direction

% global variables
epsilon = 8.854e-12; % permittivity of free space
lambda = 1/(4*pi*epsilon);
e = 1.6e-19; % charge of an electron in C
q=0.318*e; % partial charge of hydrogen
theta=52.5*pi/180; % half bond angle of H2O molecule
alpha=2*cos(theta); % dipole term
mu=-2*sin(theta)^2; %quadrupole term

pre=q*lambda;

x=[-1:0.01:5]

phi=pre.*0.5.*(1./sqrt((x+0.5*cos(52.5)).^2+sin(52.5)^2))-(1./(x-0.5*cos(52.5)))
phi_m=pre.*((2*cos(52.5))*x.^-2+0.5*(-2*sin(52.5)^2)*x.^-3)

% plot
figure(1)
hold on
plot(z,phi,'r-')
plot(z,phi_m,'b-')
legend('Exact','Multipole')