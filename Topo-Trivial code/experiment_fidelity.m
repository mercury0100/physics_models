clear;
close all;

% Import experimental data table
T = readtable('~/Dropbox/Directory/Quantum/DATA/Workbook2F.xlsx');
T = T(1:5,2:6);
Texp1 = T{:,:};

T = readtable('~/Dropbox/Directory/Quantum/DATA/Workbook5F.xlsx');
T = T(1:5,2:6);
Texp2 = T{:,:};

T = readtable('~/Dropbox/Directory/Quantum/DATA/Workbook8F.xlsx');
T = T(1:5,2:6);
Texp3 = T{:,:};

T = readtable('~/Dropbox/Directory/Quantum/DATA/Workbook3F.xlsx');
T = T(1:5,2:6);
Texp4 = T{:,:}.';

% Enter parameters
n2 = 3e-18;     % Nonlinear refractive index m^2/W
wavelength = struct('pump', 1550e-9, 'signal', 1545e-9, 'idler', 1555e-9);

w = 450e-9;     % Width of waveguide (m)
h = 220e-9;     % Height of waveguide (m)
L = 361e-6;     % Length of waveguide (m)
Aeff = w*h;     % 'Effective' area (m^2)

n = 2;          % Number of correlations to view either side of centre (even)
N = 203;        % Total number of waveguides (odd!)
ss = 1; % 1 for short-short geometry, else 0
disorder = 0;   % Disorder in the gap between waveguides
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
Tsim1 = abs(B(yvals, yvals, end)).^2;

L = 381e-6;     % Length of waveguide (m)
disorder = 0;

% Set-up Hamiltonian
[v, H, disarray, seed] = single_defect(N, ss, w, gap, wavelength, disorder/2, diag_disorder);

%% Simulate the propagation of the pump, signal and idler
P_in = 1;       % Input power
nz = 1000;      % Number of steps to compute
[A, B, C] = propagate(L, g, H, nz, P_in, ceil(N/2));

Tsim2 = abs(B(yvals, yvals, end)).^2;

L = 381e-6;     % Length of waveguide (m)
disorder = 0;

% Set-up Hamiltonian
[v, H, disarray, seed] = single_defect(N, ss, w, gap, wavelength, disorder/2, diag_disorder);

%% Simulate the propagation of the pump, signal and idler
P_in = 1;       % Input power
nz = 1000;      % Number of steps to compute
[A, B, C] = propagate(L, g, H, nz, P_in, ceil(N/2));

Tsim3 = abs(B(yvals, yvals, end)).^2;

L = 381e-6;     % Length of waveguide (m)
disorder = 0;

% Set-up Hamiltonian
[v, H, disarray, seed] = single_defect(N, ss, w, gap, wavelength, disorder/2, diag_disorder);

%% Simulate the propagation of the pump, signal and idler
P_in = 1;       % Input power
nz = 1000;      % Number of steps to compute
[A, B, C] = propagate(L, g, H, nz, P_in, ceil(N/2));

Tsim2 = abs(B(yvals, yvals, end)).^2;

L = 364e-6;     % Length of waveguide (m)
disorder = 0;

% Set-up Hamiltonian
[v, H, disarray, seed] = single_defect(N, ss, w, gap, wavelength, disorder/2, diag_disorder);

%% Simulate the propagation of the pump, signal and idler
P_in = 1;       % Input power
nz = 1000;      % Number of steps to compute
[A, B, C] = propagate(L, g, H, nz, P_in, ceil(N/2));

Tsim4 = abs(B(yvals, yvals, end)).^2;

Texp1 = Texp1./sum(Texp1,'all');
Texp2 = Texp2./sum(Texp2,'all');
Texp3 = Texp3./sum(Texp3,'all');
Texp4 = Texp4./sum(Texp4,'all');


Tsim1 = Tsim1./sum(Tsim1,'all');
Tsim2 = Tsim2./sum(Tsim2,'all');
Tsim3 = Tsim3./sum(Tsim3,'all');
Tsim4 = Tsim4./sum(Tsim4,'all');

f1 = overlap(Tsim1,Texp1)
f2 = overlap(Tsim2,Texp2);
f3 = overlap(Tsim3,Texp3);
f4 = overlap(Tsim4,Texp4);


(f1+f2+f3+f4)/4