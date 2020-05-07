function [cc, dbeta] = coupling_constant_nonidentical(gap, width1, width2, wavelength)

    % Interpolates coupling constant for given widths (m) and gap (m) 
    % arrays at an input wavelength (m). Currently, only widths of 450 nm 
    % and 465 nm are supported. Outputs the coupling constant and 
    % difference in propagation constant.
    
    dir = 'n_eff/';
    filename0 = [num2str(wavelength*10^9) '_propconst.mat'];
    filename1 = [num2str(wavelength*10^9) 'n_nonidentical_v1'];
    
    s0 = load([dir, filename0]);
    s1 = load([dir, filename1]);
    
    dn0 = interp1(s0.width_array, s0.n, width2, 'spline') - ...
        interp1(s0.width_array, s0.n, width1, 'spline');
    dn1 = interp1(s1.gap_array, s1.n1, gap, 'spline') - ...
        interp1(s1.gap_array, s1.n2, gap, 'spline');
    
    cc = pi*sqrt(dn1.^2 - dn0.^2)./s0.wavelength;
    dbeta = interp1(s0.width_array, s0.beta, width2, 'spline') - ...
        interp1(s0.width_array, s0.beta, width1, 'spline');