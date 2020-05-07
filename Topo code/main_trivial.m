clear;
close all;

% for trial = 1:8
    
% close all;

% Enter parameters
n2 = 3e-18;     % Nonlinear refractive index m^2/W
wavelength = struct('pump', 1550e-9, 'signal', 1545e-9, 'idler', 1555e-9);

w1 = 450e-9;    % Width of waveguide (m)
w2 = 465e-9;    % Width of trivial defect waveguide (m)
h = 220e-9;     % Height of waveguide (m)
L = 381e-6;     % Length of waveguide (m)
Aeff = w1*h;    % 'Effective' area (m^2)

N = 204;        % Total number of waveguides (even)
n = 10;         % Number of waveguides between trivial defects
disorder = 0.2;   % Disorder in the gap between waveguides
diag_disorder = 0e-9;   % Disorder in width of waveguides (m)

gap = 100e-9;   % Width of gap (m)

temp = 2*pi*n2/Aeff;
g = struct('pump', temp/wavelength.pump,'signal', ...
    temp/wavelength.signal, 'idler', temp/wavelength.idler);
clear('temp');

% Set-up Hamiltonian
[H, disarray, seed] = two_triv_defect(N, n, w1, w2, gap, wavelength, disorder/2, diag_disorder);

%% Simulate the propagation of the pump, signal and idler
P_in = 1;       % Input power
nz = 1000;      % Number of steps to compute
[A, B, C] = propagate(L, g, H, nz, P_in, [N/2 - n/2, N/2 + n/2 + 1]);
% [97, 108]);

%% Plotting and saving
x = linspace(1, L*10^6, nz);
y = -(N-1)/2:(N-1)/2;

dir = 'Plots/Simulations/';
folder = ['Trivial Disorder ', num2str(disorder*100), ' + ', num2str(diag_disorder*10^9), 'nm/'];
% trial = num2str(n);
% trial = num2str(input('Enter trial number: '));
% prefix = [dir, folder, num2str(trial), '_2trivial', num2str(disorder*100), '_'];
% prefix = [dir, folder, '2trivial', num2str(disorder*100), '_'];

figure(1); pcolor(x, y, abs(A(:,:)).^2); 
colormap jet; shading interp; c = colorbar; c.Label.String = 'Power (normalised)';
xlabel(['Propagation length (' char(956) 'm)']); ylabel('Waveguide number'); 
title(['Pump Propagation: ', num2str(disorder*100), '% disorder']);
axis([-inf, inf, -25, 25]);
% saveas(1, [prefix, 'pump'], 'png');

figure(2); pcolor(x, y, abs(C(:,:)).^2); 
colormap jet; shading interp; c = colorbar; c.Label.String = 'Power (normalised)';
xlabel(['Propagation length (' char(956) 'm)']); ylabel('Waveguide number'); 
title(['Signal Propagation: ', num2str(disorder*100), '% disorder']);
axis([-inf, inf, -25, 25]);
% saveas(2, [prefix, 'signal'], 'png');

% figure(3); surf(-(N-1)/2-0.5:(N-1)/2-0.5, -(N-1)/2-0.5:(N-1)/2-0.5, abs(B(:,:,end)).^2);
% axis([-105, 105, -105, 105, -inf, inf]);
% xlabel('Signal'); ylabel('Idler'); zlabel('Power (normalised)');
% view([45, 30]); title([num2str(disorder*100), '% disorder']);
% % saveas(3, [prefix, 'surfplot'], 'fig');

%%
figure(4);
% yvals = -10:10;
yvals = N/2-n:N/2+n;
% yvals = -floor(length(N/2-n:N/2+n)/2):floor(length(N/2-n:N/2+n)/2);
% b = bar3(yvals, abs(B(N/2-n:N/2+n, N/2-n:N/2+n, end)).^2);
b = bar3(yvals, abs(B(yvals, yvals, end)).^2);
for k = 1:length(b)
    b(k).CData = b(k).ZData;
    b(k).XData = b(k).XData + (min(yvals)-1)*ones(size(b(k).XData));
    b(k).Parent.YDir = 'normal';
end
colormap jet; shading interp; colorbar; axis tight;
view([45, 30]); %title([num2str(disorder*100), '% disorder']);
xlabel('Signal'); ylabel('Idler'); zlabel('Power (normalised)');
% saveas(4, [prefix, 'barplot'], 'fig');
% saveas(4, [prefix, 'barplot'], 'png');

figure(5); p1 = plot(y, abs(A(:, end).^2)); grid on;
xlabel('Waveguide Number'); ylabel('Power (normalised)'); 
title(['Output Pump Power: ', num2str(disorder*100), '% disorder']); 
p1.Parent.XLim = [-40, 40];
% saveas(5, [prefix, 'pump_output'], 'png');

figure(6); p2 = plot(y, abs(C(:, end).^2)); grid on;
xlabel('Waveguide Number'); ylabel('Power (normalised)'); 
title(['Output Signal Power: ', num2str(disorder*100), '% disorder']);
p2.Parent.XLim = [-20, 20];
% saveas(6, [prefix, 'signal_output'], 'png');

[eigvec, eigval] = eig(H.pump);
E = sum(eigval);

figure(7); plot(y, E, '.'); grid on;
hold on; plot(y(end-1:end), E(end-1:end), 'r.', 'MarkerSize', 12);
xlabel('Mode number'); ylabel('Relative $k_z$ (m$^{-1}$)'); 
% saveas(7, [prefix, 'energygap'], 'png');

%% ---------------------- HOM INTERFERENCE AT THE OUTPUT -------------------
BS = eye(N);
for k=1:N-1
    BS(N-k+1, k) = BS(N-k+1, k)+1i;
end
%TM(end,1) = 1i;
BS = BS/sqrt(2);
D = BS*B(:, :, end)*BS;

d = [B((N-n)/2, (N-n)/2, end), B((N-n)/2, (N+n)/2+1, end); ...
    B((N+n)/2+1, (N-n)/2, end), B((N+n)/2+1, (N+n)/2+1, end)];
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
b = bar3(yvals, abs(D(N/2-n:N/2+n, N/2-n:N/2+n, end)).^2);
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

peak1 = abs(D((N-n)/2, (N+n)/2 + 1).^2);
peak2 = abs(D((N+n)/2 + 1, (N-n)/2).^2);
peak3 = abs(D((N-n)/2, (N-n)/2).^2);
peak4 = abs(D((N+n)/2 + 1, (N+n)/2 + 1).^2);

vis = (peak1 + peak2 - peak3 - peak4)/(peak1 + peak2 + peak3 + peak4);

