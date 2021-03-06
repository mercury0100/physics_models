% Choose kT, N (linear lattice size), and J (strength of the coupling)
tic
kT = 3.5;
N = 30;
J = 1; % change sign for antiferromagnetic coupling

h = -1:0.1:1;
iterations =6;

Mstore = zeros(iterations,length(h));

for iter = 1:iterations
for h_n = 1:length(h)

% choose t Metropolis update steps (should be multiple of N^2)
t = 300*N^2;

% Generate a random initial configuration. 
% Comment out to keep sampling with the previous configuration.
p=.5; % average proportion of initial +1 spins
grid = sign(p-rand(N)); %random

% Run the Metropolis algorithm and return a matrix of spin values.
grid = metropolis(N,kT,J,t,grid,h(h_n));

% Compute final magnetization density and energy density
M = sum(sum(grid))/numel(grid);
E = isingenergy(grid,J);

% Plot correlation function
% cor = correlation(grid); figure(2); surf(cor); 

Mstore(iter, h_n) = M;

end
end
av_Mstore = mean(Mstore);

figure(1)
plot(h, av_Mstore);
title('h vs Magnetization for kT =' + kT);
xlabel('h');
ylabel('Average Magnetization');

toc