clear all; close all;
deb=3.33e-30;
d=1e-10;
theta=105;
thalf=theta/2;
qh=(1.86*deb)/(2*d*cosd(thalf))
qo=-2*qh
epsilon=8.85e-12;
dipole=[0 0 -2*qh*d*cosd(thalf)];
rx=0;
ry=0;
OO=[0,0,(d/2)*cosd(thalf)];
OH1=[-d*sin(thalf),0,-(d/2)*cosd(thalf)];
OH2=[d*sin(thalf),0,-(d/2)*cosd(thalf)];
quad_diag=2*d^2*sind(thalf)^2*qh*[2 , -1 , -1];
range=(0.1:0.005:0.4)*1e-9;
count=1;
for rz=range
    r=[rx,ry,rz];
    exact(count)=(4*pi*epsilon)^-1*(qo/norm(r-OO)+qh/norm(r-OH1)+qh/norm(r-OH2));
    approx(count)=(4*pi*epsilon)^-1*(dot(dipole,r)/(norm(r)^3)+0.5*dot(quad_diag,r.*r)/(norm(r)^5));
    perc(count)=abs(exact(count)-approx(count))/abs(exact(count))*100;
    count=count+1;
end

figure(1)
plot(range,perc)
title('Percent Error')

figure(2)
hold on
plot(range,exact, '-r')
plot(range,approx, '-b')
legend('Exact','Multipole')