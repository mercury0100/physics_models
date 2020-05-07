clear;
close all;

% Import experimental data table
T = readtable('~/Dropbox/Directory/Quantum/DATA/Workbook5F.xlsx');
T = T(1:5,2:6);
T = T{:,:};%.';

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
disorder = 0.2;   % Disorder in the gap between waveguides
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
yvals = (N+1)/2-n:(N+1)/2+n;

%%
figure('position', [0, 0, 1500, 1000]);

subplot(1,2,1)

b = bar3(yvals, abs(B(yvals, yvals, end)).^2);

for k = 1:length(b)
    b(k).CData = b(k).ZData;
    b(k).XData = b(k).XData + (min(yvals)-1)*ones(size(b(k).XData));
    b(k).Parent.YDir = 'normal';
end

colormap jet; shading interp; axis tight;
view([45, 30]); title(['Simulation']);
xlabel('Signal'); ylabel('Idler'); zlabel('Power (normalised)');
% saveas(4, [prefix, 'barplot'], 'fig');
% saveas(4, [prefix, 'barplot'], 'png');

subplot(1,2,2)

b = bar3(T);

for k = 1:length(b)
    b(k).CData = b(k).ZData;
    b(k).XData = b(k).XData + (min(yvals)-1)*ones(size(b(k).XData));
    b(k).Parent.YDir = 'normal';
end

colormap jet; shading interp; axis tight;
view([45, 30]); title(['Experiment']);
yticks([1 2 3 4 5]);
yticklabels({'100','101','102','103', '104'});
xlabel('Signal'); ylabel('Idler'); zlabel('Counts');

