function [cc, lc] = coupling_constant(width, gap, wavelength)

    % Interpolates coupling length for given width (m) and gap (m) arrays
    % at an input wavelength (m). Outputs the coupling constants and
    % coupling lengths in metres in an array where entry ij is the ith 
    % width and jth gap.
    
    dir = 'n_eff/';
    filename = [num2str(wavelength*10^9), 'dataset_v2.mat'];
    s = load([dir, filename]);
    
    [Xq, Yq] = meshgrid(width, gap);
    lc = interp2(s.width_array, s.gap_array, real(s.lc), Xq, Yq, 'spline');
    
    cc = pi./(2*lc);