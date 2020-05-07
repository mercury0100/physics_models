function [v, H, disarray, seed] = two_topo_defect(N, n, width, gap, wavelength, disorder, diag_disorder)

    % Calculates the Hamiltonian of the system. Takes as input:
    %    N = total number of waveguides
    %    n = number of waveguides between the two topological defects
    %    width = width of waveguides
    %    gap = gap(s) between the waveguides (without disorder)
    %    disorder = amount of disorder to add
    %    and gives as output:
    %    H = Hamiltonian for the pump, signal and idler as a structure

    sep = n/2;
    u = zeros(1, sep);
    v = zeros(1, N/2 - 1);
    
    gap_long = max(gap);
    gap_short = min(gap);

    for k = 1:N/2-sep-1
        if mod(k, 2) == 1
            v(k+sep) = gap_long;
        else
            v(k+sep) = gap_short;
        end
    end
    
    for k = 1:sep
        if mod(k, 2) == 1
            u(k) = gap_long;
        else
            u(k) = gap_short;
        end
    end

    v(1:sep) = flip(u);

    if mod(sep, 2) == 1
        v = [flip(v), gap_short, v];
    else
        v = [flip(v), gap_long, v];
    end
    
    cc.pump = coupling_constant(width, v, wavelength.pump);
    cc.signal = coupling_constant(width, v, wavelength.signal);
    cc.idler = coupling_constant(width, v, wavelength.idler);
    
    % Add disorder to coupling constants (waveguide gaps)
    seed = rng;
    if disorder ~= 0
        %disarray = disorder - 2*disorder*rand(length(v), 1);
        disarray = disorder*normrnd(0, disorder, length(v), 1);
        cc.pump = cc.pump + disarray.*cc.pump;
        cc.signal = cc.signal + disarray.*cc.signal;
        cc.idler = cc.idler + disarray.*cc.idler;
    else
        disarray = 0;
    end
    
    % Add disorder to propagation constants (waveguide widths)
    dbeta = calc_dbeta(N, width, diag_disorder);
        
    % Create Hamiltonian
    H.pump = diag(cc.pump, 1) + diag(dbeta.pump) + diag(cc.pump, -1);
    H.signal = diag(cc.signal, 1) + diag(dbeta.signal) + diag(cc.signal, -1);
    H.idler = diag(cc.idler, 1) + diag(dbeta.idler) + diag(cc.idler, -1);
    
    