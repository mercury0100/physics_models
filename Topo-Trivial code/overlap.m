function [OI] = overlap(x,y)
    M=5;
    xc = linspace(1,5,M);
    yc = linspace(1,5,M);
    A = x .* y;
    B = trapz(xc,trapz(yc, A, 1))^2;
    C = trapz(xc,trapz(yc, x.^2, 1)) * trapz(xc,trapz(yc, y.^2, 1));
    OI = B./C;
    OI
end