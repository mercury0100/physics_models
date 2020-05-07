% Recursive algorithm generating a tree to logical solutions of the Collatz
% conjecture based on the trivial solution, n=1. i.e. For any number x that
% satisfies the conjecture, 2x also satisfies the conjecture. If x is even,
% then (x-1)/3 is also a solution, if it exists.
clear all
format short

x(1)=8;
i=1;
it=20;
list=zeros(1);

G=graph

while i <= it
    l=[];
    for n=1:length(x)
        if rem(x(n),2)==0
            if rem((x(n)-1)/3,1)==0
                y=[2*x(n),(x(n)-1)/3];
            else
                y=[2*x(n)];
            end
        else
            y=[2*x(n)];
        end
        l=[l,y];
    end
    list(i,1:length(l))=l;
    x=l;
    i=i+1;
end
        
        