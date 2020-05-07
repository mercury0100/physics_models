a=[1 1 1];
b=[1 1 1];
n=conv(a,b);
d=[1 0 1 1];
[q,r]=deconv(n,d);
r