function [H, disarray, seed] = two_triv_defect(N, n, width1, width2, gap, wavelength, disorder, diag_disorder)

    % Calculates the Hamiltonian of the system. Takes as input:
    %    N = total number of waveguides
    %    n = number of waveguides between the two trivial defects
    %    width1 = width of waveguides
    %    width2 = width of the defect (greater than width1)
    %    gap = gap(s) between the waveguides (without disorder)
    %    disorder = amount of disorder to add
    % and gives as output:
    %    H = Hamiltonian for the pump, signal and idler as a structure
    
    c0.pump = coupling_constant(width1, gap, wavelength.pump);
    c0.signal = coupling_constant(width1, gap, wavelength.signal);
    c0.idler = coupling_constant(width1, gap, wavelength.idler);

    [c1.pump, dbeta.pump] = coupling_constant_nonidentical(gap, width1, width2, wavelength.pump);
    [c1.signal, dbeta.signal] = coupling_constant_nonidentical(gap, width1, width2, wavelength.signal);
    [c1.idler, dbeta.idler] = coupling_constant_nonidentical(gap, width1, width2, wavelength.idler);

    cc.pump = c0.pump*ones(1, N/2 - 1);
    cc.pump(n/2) = c1.pump; cc.pump(n/2 + 1) = c1.pump;
    cc.pump = [flip(cc.pump), c0.pump, cc.pump];
    
    cc.signal = c0.signal*ones(1, N/2 - 1);
    cc.signal(n/2) = c1.signal; cc.signal(n/2 + 1) = c1.signal;
    cc.signal = [flip(cc.signal), c0.signal, cc.signal];
    
    cc.idler = c0.idler*ones(1, N/2 - 1);
    cc.idler(n/2) = c1.idler; cc.idler(n/2 + 1) = c1.idler;
    cc.idler = [flip(cc.idler), c0.idler, cc.idler];
    
    % Add disorder to coupling constants
    seed = rng;
    if disorder ~= 0
        disarray = disorder - 2*disorder*rand(1, length(cc.pump));
        cc.pump = cc.pump + disarray.*cc.pump;
        cc.signal = cc.signal + disarray.*cc.signal;
        cc.idler = cc.idler + disarray.*cc.idler;
    else
        disarray = 0;
    end
    
    % Add disorder to propagation constants
    db = calc_dbeta(N, width1, diag_disorder);
    
    H.pump = diag(cc.pump, 1) + diag(db.pump) + diag(cc.pump, -1);
    H.pump(N/2 - n/2, N/2 - n/2) = H.pump(N/2 - n/2, N/2 - n/2) + dbeta.pump;
    H.pump(N/2 + n/2 + 1, N/2 + n/2 + 1) = H.pump(N/2 + n/2 + 1, N/2 + n/2 + 1) + dbeta.pump;
    
    H.signal = diag(cc.signal, 1) + diag(db.signal) + diag(cc.signal, -1);
    H.signal(N/2 - n/2, N/2 - n/2) = H.signal(N/2 - n/2, N/2 - n/2) + dbeta.signal;
    H.signal(N/2 + n/2 + 1, N/2 + n/2 + 1) = H.signal(N/2 + n/2 + 1, N/2 + n/2 + 1) + dbeta.signal;
    
    H.idler = diag(cc.idler, 1) + diag(db.idler) + diag(cc.idler, -1);
    H.idler(N/2 - n/2, N/2 - n/2) = H.idler(N/2 - n/2, N/2 - n/2) + dbeta.idler;
    H.idler(N/2 + n/2 + 1, N/2 + n/2 + 1) = H.idler(N/2 + n/2 + 1, N/2 + n/2 + 1) + dbeta.idler;
    