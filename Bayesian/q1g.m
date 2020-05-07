close all; 
clear all;
format short

%define mean, data and physical constants
vo = 1.65;
g = 9.8*10^-3;
Data = [51,198,147,39,144,288,226,13,89,137,98,106,102,100,55,44,129,22,77,82];

OddsRatio = 1;
E1=0;
E2=0;

dx=0.1;
dsig=0.01;
X=[50:dx:250];
sig=[0.05:dsig:0.6];

%integrate over all of parameter space
for x=1:length(X)
    for y=1:length(sig)
        PriA = 1;
        PriB = 1;
        for i=1:length(Data)
            ri = Data(i);
            PriA = A(PriA,sig(y),X(x),ri,g,vo);
            PriB = B(PriB,sig(y),X(x),ri,g,vo);
        end
        E1 = E1 + PriA;
        E2 = E2 + PriB;
        Odds = E2/E1
    end
end

ProbabilityOfModel2 = 1/(1+1/Odds)

function likeli = A(Prior,Sig,X,ri,g,vo)
    vi = sqrt(g.*(ri+X));
    UpperExp = (-(vi - vo).^2)./(2*Sig.^2);
    BottomFrac = 2*Sig.*sqrt(2*pi*g*(ri+X));
    like = (g.*exp(UpperExp))./BottomFrac;
    likeli= Prior.*like;
end

function likelif = B(Prior,Sig,X,ri,g,vo)
    vi = g*(ri+X)/1.9;
    UpperExp = (-(vi - vo).^2)./(2*Sig.^2);
    BottomFrac = 1.9*sqrt(2*pi)*Sig;
    like = (g.*exp(UpperExp))./BottomFrac;
    likelif= Prior.*like;
end