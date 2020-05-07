clear;
close all;

%% Set up the tight binding Hamiltonian
load('n_eff/silicon_dimer_var.mat');
% load('silicon_dimer_var_500um.mat');
plot(gap, coupling_length); grid on;
xlabel('Gap (m)'); ylabel('Coupling Constant (\mum)');
title('Silicon Waveguides 450 nm by 220 nm');
% title('Silicon Waveguides 500 nm by 220 nm');

defect = 'Short';
N = 203;
core_width = 450e-9;
v = zeros(1, N-1);

% Input gap distances
if strcmp(defect, 'Long')
    % Long-long defect --> da' = db = 166 nm and da = db' = 294 nm
    gap_short = 166e-9;
    gap_long = 294e-9;
elseif strcmp(defect, 'Equal')
    % Equidistant waveguides
    gap_short = 230e-9;
    gap_long = 230e-9;  
elseif strcmp(defect, 'Short') 
    % Short-short defect --> da' = db = 182 nm and da = db' = 324 nm
    gap_short = 182e-9;
    gap_long = 324e-9;
end

% Calculate coupling length and coupling constant
[p, s, mu] = polyfit(gap, coupling_length, 10);
cl_short = polyval(p, gap_short, s, mu);
cl_long = polyval(p, gap_long, s, mu);

c_short = pi/(2*cl_short);
c_long = pi/(2*cl_long);
% hold on; plot(gap_short, cl_short, 'rx', gap_long, cl_long, 'rx');

% Set up the Hamiltonian matrix
if strcmp(defect, 'Long') || strcmp(defect, 'Equal') 
    for k = 1:length(v)
        if k < 0.5*length(v)
            if mod(k, 2) == 1
                v(k) = c_short;
            else
                v(k) = c_long;
            end
        else
            if mod(k, 2) == 1
                v(k) = c_long;
            else
                v(k) = c_short;
            end
        end
    end
    
elseif strcmp(defect, 'Short')
    for k = 1:length(v)
        if k < 0.5*length(v) + 1
            if mod(k, 2) == 1
                v(k) = c_short;
            else
                v(k) = c_long;
            end
        else
            if mod(k, 2) == 1
                v(k) = c_long;
            else
                v(k) = c_short;
            end
        end
    end
    
end

H = 2*(diag(v, 1) + diag(v, -1));
[eigvec, eigval] = eig(H);
E = sum(eigval);

if strcmp(defect, 'Equal')
    fprintf('\nEquidistant waveguides: \n');
else
    fprintf('\n%s-%s defect: \n', defect, defect);
end
fprintf('t_a''/t_a = %.4f \n', c_short/c_long);
fprintf('t_b''/t_b = %.4f \n\n', c_long/c_short);

%% Simulate the propagation of the pump, signal and idler

%**************************************************************************
%----------------------- FIXED PARAMETERS ---------------------------------
% Inputs
% s = N;
% n1 = (s-3)/2;
L = 381e-6;     % Length of structure
% L = 500e-6;
nz = 1e3;       % Number of steps along length
dz = L/nz;      % Step size
g = 120;        % Nonlinear coefficient of WG, SPM per watt per meter (250 for significant effects)
P_in = 1;       % Input peak power
H_p = H;
H_s = H;
H_i = H;

%************************** INPUT FIELD ***********************************
% Input pump field
A = zeros(N,nz);
% A(floor(N/2), 1) = sqrt(P_in);
% A(floor(N/2) , 1) = sqrt(P_in);
A(ceil(N/2), 1) = sqrt(P_in);

% Signal/idler field
B = zeros(N,N,nz);
C = zeros(N,nz);

%**************************************************************************
%*************************** MAIN LOOP ************************************
%**************************************************************************
for j=1:nz-1
    A(:,j+1) = A(:,j) + 1i*H_p*dz*A(:,j) - 0.5*dz^2*H_p^2*A(:,j);
    B(:,:,j+1) = B(:,:,j) + 1i*H_s*dz*B(:,:,j) + 1i*dz*B(:,:,j)*H_i - 0.5*dz^2*H_s^2*B(:,:,j)...
        - 0.5*dz^2*B(:,:,j)*H_i^2-dz^2*H_s*B(:,:,j)*H_i + 1i*g*dz*diag(A(:,j).^2);
    C(:,j+1) = diag(B(:,:,j+1));
end

%% Plotting
x = linspace(1, L*10^6, nz);
y = -(N-1)/2:(N-1)/2;

figure(1); pcolor(x, y, abs(A(:,:)).^2); 
colormap jet; shading interp; c = colorbar; c.Label.String = 'Power (normalised)';
xlabel(['Propagation length (' char(956) 'm)']); ylabel('Waveguide number'); 
title(['Pump Propagation: ', defect]);
axis([-inf, inf, -40, 40]);

figure(2); pcolor(x, y, abs(C(:,:)).^2); 
colormap jet; shading interp; c = colorbar; c.Label.String = 'Power (normalised)';
xlabel(['Propagation length (' char(956) 'm)']); ylabel('Waveguide number'); 
title(['Signal Propagation: ', defect]);
axis([-inf, inf, -40, 40]);

figure(3); 
surf(-(N-1)/2+0.5:(N-1)/2+0.5, -(N-1)/2+0.5:(N-1)/2+0.5, abs(B(:,:,end)).^2);
axis([-105, 105, -105, 105, -inf, inf]);
xlabel('Signal'); ylabel('Idler'); title(defect);

figure(4);
yvals = -10:10;
% yvals = -floor(length(N/2-n:N/2+n)/2):floor(length(N/2-n:N/2+n)/2);
% b = bar3(yvals, abs(B(N/2-n:N/2+n, N/2-n:N/2+n, end)).^2);
b = bar3(yvals, abs(B(92:112, 92:112, end)).^2);
for k = 1:length(b)
    b(k).CData = b(k).ZData;
    b(k).XData = b(k).XData + (min(yvals)-1)*ones(size(b(k).XData));
    b(k).Parent.YDir = 'normal';
end
colormap jet; shading interp; colorbar; axis tight;
view([45, 30]); title(defect);
xlabel('Signal'); ylabel('Idler'); zlabel('Power (normalised)');

figure(5); p1 = bar(y, abs(A(:, end)).^2); grid on;
xlabel('Waveguide Number'); ylabel('Power (normalised)'); 
title(['Output Pump Power: ', defect]); 
p1.Parent.XLim = [-40, 40];

figure(6); p2 = bar(y, abs(C(:, end).^2)); grid on;
xlabel('Waveguide Number'); ylabel('Power (normalised)'); 
title(['Output Signal Power: ', defect]);
p2.Parent.XLim = [-20, 20];

figure(7);
plot(y, E, '.'); grid on;
xlabel('Mode number'); ylabel('k (m^{-1})');
title(['Energy Levels: ' defect]);