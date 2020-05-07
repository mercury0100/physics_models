clear;
close all;

% Enter parameters
n2 = 3e-18;     % Nonlinear refractive index m^2/W
wavelength = struct('pump', 1550e-9, 'signal', 1545e-9, 'idler', 1555e-9);

w = 450e-9;     % Width of waveguide (m)
h = 220e-9;     % Height of waveguide (m)
L = 381e-6;     % Length of waveguide (m)
Aeff = w*h;     % 'Effective' area (m^2)

n = 2;          % Number of correlations to view either side of centre (even)
N = 203;        % Total number of waveguides (odd!)
ss = 1; % 1 for short-short geometry, else 0
disorder = 0.8;   % Disorder in the gap between waveguides
diag_disorder = 0; %1e-9;   % Disorder in width of waveguides (m)

gap_short = 173e-9;             % Width of short gap (m)
gap_long = 307e-9;              % Width of long gap (m)
gap = [gap_short, gap_long];

temp = 2*pi*n2/Aeff;
g = struct('pump', temp/wavelength.pump,'signal', ...
    temp/wavelength.signal, 'idler', temp/wavelength.idler);
clear('temp');

% Set-up Hamiltonian
[v, H, disarray, seed] = single_defect(N, ss, w, gap, wavelength, disorder/2, diag_disorder);

%% Simulate the propagation of the pump, signal and idler
P_in = 1;       % Input power
nz = 1000;      % Number of steps to compute
[A, B, C] = propagate(L, g, H, nz, P_in, ceil(N/2));
% [97, 108]);

%% Plotting and saving
x = linspace(1, L*10^6, nz);
y = -(N-1)/2:(N-1)/2;

% dir = 'Plots/Simulations/';
% folder = ['Disorder ', num2str(disorder*100), ' + ', num2str(diag_disorder*10^9), 'nm/'];
% trial = num2str(input('Enter trial number: '));
% trial = num2str(n);
% prefix = [dir, folder, num2str(trial), '_2topo', num2str(disorder*100), '_'];

%% Fig 1

figure(1); pcolor(x, y, abs(A(:,:)).^2); 
colormap jet; shading interp; c = colorbar; c.Label.String = 'Power (normalised)';
xlabel(['Propagation length (' char(956) 'm)']); ylabel('Waveguide number'); 
title(['Pump Propagation: ', num2str(disorder*100), '% disorder']);
axis([-inf, inf, -n, n]);
% saveas(1, [prefix, 'pump'], 'png');
%% Fig 2

figure(2); pcolor(x, y, abs(C(:,:)).^2); 
colormap jet; shading interp; c = colorbar; c.Label.String = 'Power (normalised)';
xlabel(['Propagation length (' char(956) 'm)']); ylabel('Waveguide number'); 
title(['Photon Pair Propagation: ', num2str(disorder*100), '% disorder']);
axis([-inf, inf, -n, n]);
% saveas(2, [prefix, 'signal'], 'png');

% figure(3); surf(-(N-1)/2-0.5:(N-1)/2-0.5, -(N-1)/2-0.5:(N-1)/2-0.5, abs(B(:,:,end)).^2);
% axis([-105, 105, -105, 105, -inf, inf]);
% xlabel('Signal'); ylabel('Idler'); zlabel('Power (normalised)');
% view([45, 30]); title([num2str(disorder*100), '% disorder']);
% saveas(3, [prefix, 'surfplot'], 'fig');

%%
figure(4);
% yvals = -10:10;
% yvals = -n:n;
yvals = (N+1)/2-n:(N+1)/2+n;
% b = bar3(yvals, abs(B(N/2-n:N/2+n, N/2-n:N/2+n, end)).^2);
% b = bar3(yvals, abs(B(92:112, 92:112, end)).^2);
b = bar3(yvals, abs(B(yvals, yvals, end)).^2);
for k = 1:length(b)
    b(k).CData = b(k).ZData;
    b(k).XData = b(k).XData + (min(yvals)-1)*ones(size(b(k).XData));
    b(k).Parent.YDir = 'normal';
end
colormap jet; shading interp; colorbar; axis tight;
view([45, 30]); title(['Disorder: ', num2str(disorder*100), '%']);
xlabel('Signal'); ylabel('Idler'); zlabel('Power (normalised)');

writematrix(abs(B(yvals, yvals, end)).^2, ['sim_output_', num2str(disorder), '.csv'])

% saveas(4, [prefix, 'barplot'], 'fig');
% saveas(4, [prefix, 'barplot'], 'png');

%%

figure(5); p1 = bar(y, abs(A(:, end).^2)); grid on;
xlabel('Waveguide Number'); ylabel('Power (normalised)'); 
title(['Output Pump Power: ', num2str(disorder*100), '% disorder']); 
p1.Parent.XLim = [-40, 40];
% saveas(5, [prefix, 'pump_output'], 'png');

%%

figure(6); p2 = bar(y, abs(C(:, end).^2)); grid on;
xlabel('Waveguide Number'); ylabel('Power (normalised)'); 
title(['Output Signal Power: ', num2str(disorder*100), '% disorder']);
p2.Parent.XLim = [-40, 40];
% saveas(6, [prefix, 'signal_output'], 'png');

[eigvec, eigval] = eig(H.pump);
E = sum(eigval);

%%

figure(7); plot(y, E, '.'); grid on;
hold on; plot(y((N+1)/2), E((N+1)/2), 'r.', 'MarkerSize', 12);
xlabel('Mode number'); ylabel({'Relative $k_z$ (m$^{-1}$)'},'Interpreter','latex'); 
title(['Disorder: ', num2str(disorder*100), '%']);
% saveas(7, [prefix, 'energygap'], 'png');

%{
%% ---------------------- HOM INTERFERENCE AT THE OUTPUT -------------------
BS = eye(N);
for k=1:N-1
    BS(N-k+1, k) = BS(N-k+1, k) + 1i;
end
BS(1, end) = 1i;
BS = BS/sqrt(2);

D = BS*B(:, :, end)*BS;

d = [B((N+1-n)/2, (N+1-n)/2, end), B((N+1-n)/2, (N+1+n)/2+1, end); ...
    B((N+1+n)/2+1, (N+1-n)/2, end), B((N+1+n)/2+1, (N+1+n)/2+1, end)];
BS2 = [1, 1i; 1i, 1];
theta = -pi:pi/128:pi;
iv = zeros(1, length(k));
for k = 1:length(theta)
    PM = [1, 0; 0, exp(1i*theta(k))];
    rho = BS2*PM*d*PM*BS2;
    difference = abs(rho(1, 2).^2) + abs(rho(2, 1).^2) - abs(rho(1, 1).^2) - abs(rho(2, 2).^2);
    sum = abs(rho(1, 2).^2) + abs(rho(2, 1).^2) + abs(rho(1, 1).^2) + abs(rho(2, 2).^2);
    iv(k) = difference/sum;
end
% figure(8); plot(theta, iv);

figure(9);
b = bar3(yvals, abs(D(yvals, yvals, end)).^2);
for k = 1:length(b)
    b(k).CData = b(k).ZData;
    b(k).XData = b(k).XData + (min(yvals)-1)*ones(size(b(k).XData));
    b(k).Parent.YDir = 'normal';
end
colormap jet; shading interp; colorbar; axis tight;
view([45, 30]); title([num2str(disorder*100), '% disorder']);
xlabel('Signal'); ylabel('Idler'); zlabel('Power (normalised)');
% saveas(9, [prefix, 'barplot'], 'fig');
% saveas(9, [prefix, 'barplot'], 'png');

peak1 = abs(D((N+1-n)/2, (N+1+n)/2 + 1).^2);
peak2 = abs(D((N+1+n)/2 + 1, (N+1-n)/2).^2);
peak3 = abs(D((N+1-n)/2, (N+1-n)/2).^2);
peak4 = abs(D((N+1+n)/2 + 1, (N+1+n)/2 + 1).^2);

vis = (peak1 + peak2 - peak3 - peak4)/(peak1 + peak2 + peak3 + peak4);
difference = abs(eigval(102, 102) - eigval(103, 103));
%}


