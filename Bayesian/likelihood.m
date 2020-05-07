
function L = likelihood(Prior,Sig,X,ri,g,vo)
    vi = sqrt(g.*(ri+X));
    Upper = g.*exp((-(vi - vo).^2)./(2*Sig.^2));
    Lower = 2*Sig.*sqrt(2*pi*g*(ri+X));
    llh = Upper./Lower;
    L = Prior.*llh;
end