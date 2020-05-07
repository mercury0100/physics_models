%Likelihood function for the Bernoulli distribution

function L = likelihood_B(Prior,Sig,X,ri,g,vo)
    vi = g*(ri+X)/1.9;
    Upper = g.*exp((-(vi - vo).^2)./(2*Sig.^2));
    Lower = 2*Sig.*sqrt(2*pi*g*(ri+X));
    llh = Upper./Lower;
    L = Prior.*llh;
end