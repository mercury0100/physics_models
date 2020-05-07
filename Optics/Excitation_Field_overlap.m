clear all
format short

%Global variables & step size
d=4e-6;
y=[-2:0.01:2]*10^-6;

%Excitation field s(y)
s=sqrt(2/d).*exp(-(y.^2)./(d/8)^2);

%Transverse field distribution u
for m = 1:5
    if rem(m,2)==0
        u(m,:)=sqrt(2/d).*sin(m*pi*y./d);
    else
        u(m,:)=sqrt(2/d).*cos(m*pi*y./d);
    end
end

%Integrate
for n = 1:5
    a(n)=dot(u(n,:),s)*0.01e-6;
end

a

figure(1)
subplot(3,2,1)
plot(u(1,:),y,s,y)
title('m=1')
subplot(3,2,2)
plot(u(2,:),y,s,y)
title('m=2')
subplot(3,2,3)
plot(u(3,:),y,s,y)
title('m=3')
subplot(3,2,4)
plot(u(4,:),y,s,y)
title('m=4')
subplot(3,2,5)
plot(u(5,:),y,s,y)
title('m=5')
subplot(3,2,6)
plot(u(1,:),y,u(2,:),y,u(3,:),y,u(4,:),y,u(5,:),y)
title('all modes')