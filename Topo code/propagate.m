function [A, B, C] = propagate(L, g, H, nz, P_in, varargin)

    % Computes propagation of input though waveguide array. Takes as input:
    %   L = length of waveguides
    %   g = nonlinear parameter
    %   H = structure containing the pump, signal and idler Hamiltonians
    %   nz = number of steps to compute
    %   P_in = input power
    %   v = vector containing waveguide number(s) of light input - default
    %   is one input at ceil(N/2)
    % and returns three variables:
    %   A = N by nz array of pump power per waveguide per time step
    %   B = N by N by nz array of signal/idler correlation per time step
    %   C = N by nz array of signal/idler power per waveguide per time step
    
    % Input variables - step size and Hamiltonians
    dz = L/nz;
    H_p = H.pump;
    H_s = H.signal;
    H_i = H.idler;
    N = length(H.pump);
    
    if nargin == 5
        v = ceil(N/2);
    else
        v = varargin{1};
    end
    
    % Input pump field
    A = zeros(N, nz);
    A(v, 1) = sqrt(P_in/length(v));
    
    % Signal/idler field
    B = zeros(N, N, nz);
    C = zeros(N, nz);
    
    % Main loop
    for j=1:nz-1
        A(:,j+1) = A(:,j) + 1i*H_p*dz*A(:,j) - 0.5*dz^2*H_p^2*A(:,j);
        B(:,:,j+1) = B(:,:,j) + 1i*H_s*dz*B(:,:,j)...
            + 1i*dz*B(:,:,j)*H_i - 0.5*dz^2*H_s^2*B(:,:,j)...
            - 0.5*dz^2*B(:,:,j)*H_i^2 - dz^2*H_s*B(:,:,j)*H_i...
            + 1i*g.pump*dz*diag(A(:,j).^2);
        C(:,j+1) = diag(B(:,:,j+1));
    end
