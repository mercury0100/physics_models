%Joint posterior space of the Bernoulli model given the data

function Z = pdf_B(X,Y)
    vo = 1.65;
    g = 9.8*10^-3;
    Data = [51,198,147,39,144,288,226,13,89,137,98,106,102,100,55,44,129,22,77,82];
    Prior = ones;
    for i=1:length(Data)
        ri = Data(i);
        Posterior=likelihood_B(Prior,Y,X,ri,g,vo);
        Prior = Posterior;
    end 
    Z = Posterior;
end