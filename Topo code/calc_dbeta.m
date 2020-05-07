function dbeta = calc_dbeta(N, width, diag_disorder)

    dw1 = -0.5*diag_disorder + diag_disorder*rand(1, N);
    dw2 = -0.5*diag_disorder + diag_disorder*rand(1, N);
    %dw1 = diag_disorder*normrnd(0,diag_disorder,1, N);
    %dw2 = diag_disorder*normrnd(0,diag_disorder,1, N);
    w = width + dw1 + dw2;
    
%     w = (width - diag_disorder) + 2*diag_disorder*rand(1, N);
%     w = normrnd(width, diag_disorder/2, 1, N);
%     w = width + diag_disorder; % Max width
%     w = width - diag_disorder; % Min width

    s_p = load('n_eff/1550_propconst.mat');
    s_s = load('n_eff/1545_propconst.mat');
    s_i = load('n_eff/1555_propconst.mat');
    
    beta.pump = interp1(s_p.width_array, s_p.beta, width, 'spline');
%     beta.signal = interp1(s_s.width_array, s_s.beta, width, 'spline');
%     beta.idler = interp1(s_i.width_array, s_i.beta, width, 'spline');

    dbeta.pump = interp1(s_p.width_array, s_p.beta, w, 'spline') - beta.pump;
    dbeta.signal = interp1(s_s.width_array, s_s.beta, w, 'spline') - beta.pump;
    dbeta.idler = interp1(s_i.width_array, s_i.beta, w, 'spline') - beta.pump;
