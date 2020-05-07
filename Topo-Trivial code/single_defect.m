function [v, H, disarray, seed] = single_defect(N, ss, width, gap, wavelength, disorder, diag_disorder)

    % Calculates the Hamiltonian of the system. Takes as input:
    %    N = total number of waveguides
    %    n = number of waveguides between the two topological defects
    %    width = width of waveguides
    %    gap = gap(s) between the waveguides (without disorder)
    %    disorder = amount of disorder to add
    %    and gives as output:
    %    H = Hamiltonian for the pump, signal and idler as a structure

    %v = zeros(1, N/2 - 1);
    v = zeros(1, (N-1)/2);
    
    gap_long = max(gap);
    gap_short = min(gap);

    if mod((N-1)/2,2)==0
        for k = 1:(N-1)/2
            if mod(k, 2) == 1
                v(k) = gap_long;
            else
                v(k) = gap_short;
            end
        end
    else
        for k = 1:(N-1)/2
            if mod(k, 2) == 1
                v(k) = gap_short;
            else
                v(k) = gap_long;
            end
        end
    end
        

    v = [v, flip(v)];
    
    cc.pump = coupling_constant(width, v, wavelength.pump);
    cc.signal = coupling_constant(width, v, wavelength.signal);
    cc.idler = coupling_constant(width, v, wavelength.idler);
    
    % Add disorder to coupling constants (waveguide gaps)
    seed = rng;
    if disorder ~= 0
        %disarray = disorder - 2*disorder*rand(length(v), 1);
        disarray = disorder*normrnd(0, 1/sqrt(12), length(v), 1);
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
    
    