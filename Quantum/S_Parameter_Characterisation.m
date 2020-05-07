% Characterisation of passive components using S parameters

%Given
R=20000;
Zo=50;
f=1e9;
S11_val=0.333;
S11_ph=30;

%Conversion to cartesian
S11_re=S11_val*cos(S11_ph*pi/180);
S11_im=S11_val*sin(S11_ph*pi/180);

%Numerator Addition
Z_num_re=1+S11_re;
Z_num_im=S11_im;

%Denominator Addition
Z_den_re=1-S11_re;
Z_den_im=S11_im*(-1);

%Conversion to polar
Z_num_val=sqrt(Z_num_re^2+Z_num_im^2);
Z_num_ph=atan(Z_num_im/Z_num_re);
Z_den_val=sqrt(Z_den_re^2+Z_den_im^2);
Z_den_ph=atan(Z_den_im/Z_den_re);

%Division
Zq_val=Z_num_val/Z_den_val;
Zq_ph=Z_num_ph-Z_den_ph;

%Multiplication by Zo
Z_val=Zo*Zq_val;
Z_ph=Zq_ph*180/pi;

%Back to Cartesian
Z_re=Z_val*cos(Z_ph);
Z_im=Z_val*sin(Z_ph);

%Capacitor Dissipation Factor and Capacitance
C=1/(2*pi*f);
D=Z_re/Z_im;

%Inductor Quality factor and inductance
L=Z_im/(2*pi*f);
Q=Z_im/Z_re;
