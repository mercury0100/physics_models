clear all
format short

%Globals
d=4e-6;
x=[-2:0.01:2]*10^-6;

for m = 1:5
    if rem(m,2)==0
        u(m,:)=sqrt(2/d).*sin(m*pi*x./d);
    else
        u(m,:)=sqrt(2/d).*cos(m*pi*x./d);
    end
end

plot(u(1,:),x,u(2,:),x,u(3,:),x,u(4,:),x,u(5,:),x)