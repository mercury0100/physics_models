close all
clear all

%---- Directory for saving the figures --------------------------------
directory = '/Users/Cooper/Desktop/TOPOLOGICAL PHOTONICS/';
SAVEFIGS =0;
%inputs
N=203; %number of waveguides. Currently needs to be 1 greater than a multiple of 4
s=N;
n1=(s-3)/2;
L=0.6E-3; %length of structure
nz = 1E3; %number of steps along length
dz=L/nz; %size of a step
g=250; %nonlinear coefficient of wg, SPM per watt per meter
P_in=1; %input peak power

%wavelengths
l_pump=1550E-9; %pump wavelength
l_signal=1545E-9; %signal wavelength
l_idler=(2/l_pump-1/l_signal)^-1; %idler wavelength for frequency matching

% ********************** CHOSING DESIGN **********************************
d_2_7 =0;
d_3_7 =1;
d_2_4 =0;
d_3_4 =0;
swap_coupling_constants=0; %if this is 1, swap over the coupling constants, for a short-short defect
% *************************************************************************

% ********************** INTRODUCING DEFECTS ******************************
disorder                = 0;        % Disorder in the position of the waveguides (off-diagonal)
noise_lvl               = 0.75e4;    % Level of disorder (off-diagonal)
percent                 = 0;

missing_waveguides      = 0;        % Missing waveguides (off-diagonal)
defect_position         = 105;      % Position of missing waveguide

wider_center_waveguide  = 0;        % Wider central waveguide in an array of equidistant waveguides (on-diagonal)--> It creates a trivial defect
w_factor                = 1.5;      % Widening factor

on_diagonal             = 0;
cc_noise_lvl            = 1.5*2.77e4; 
% *************************************************************************

%***************************Coupling constants****************************
coupling    = load('silicon_dimer_var.mat'); % Load Coupling curve
if d_2_7
    topo_quantum_sim_coupling_2_7; %load coupling constants as a function of wavelength
    C_close_p = spline(ll,C_close,l_pump); %coupling at the pump wavelength
    C_far_p = spline(ll,C_far,l_pump);
    C_close_s = spline(ll,C_close,l_signal); %coupling at the signal wavelength
    C_far_s = spline(ll,C_far,l_signal);
    C_close_i = spline(ll,C_close,l_idler); %coupling at the idler wavelength
    C_far_i = spline(ll,C_far,l_idler);
   
end
if d_3_7
    topo_quantum_sim_coupling_3_7; %load coupling constants as a function of wavelength
    %C_close = C_far;
    C_close_p = spline(ll,C_close,l_pump); %coupling at the pump wavelength
    C_far_p = spline(ll,C_far,l_pump);
    C_close_s = spline(ll,C_close,l_signal); %coupling at the signal wavelength
    C_far_s = spline(ll,C_far,l_signal);
    C_close_i = spline(ll,C_close,l_idler); %coupling at the idler wavelength
    C_far_i = spline(ll,C_far,l_idler);
end

if swap_coupling_constants==1
    x=C_close_p;
    C_close_p=C_far_p;
    C_far_p=x;
    x=C_close_s;
    C_close_s=C_far_s;
    C_far_s=x;
    x=C_close_i;
    C_close_i=C_far_i;
    C_far_i=x;
end

if d_2_4
    C_close_p = 27686.55203; %coupling at the pump wavelength
    C_far_p = C_close_p;
    C_close_s = 26891.62644; %coupling at the signal wavelength
    C_far_s = C_close_s;
    C_close_i = 28504.65; %coupling at the idler wavelength
    C_far_i = C_close_i;
end
if d_3_4
    C_close_p = 24792.23312; %coupling at the pump wavelength
    C_far_p = C_close_p;
    C_close_s = 24057.07617; %coupling at the signal wavelength
    C_far_s = C_close_s;
    C_close_i = 25548.92649; %coupling at the idler wavelength
    C_far_i = C_close_i;
end
%*************************************************************************

%***********Nearest neighbour coupling Hamiltonian************************
array1p = [];
array1s = [];
array1i = [];
for ii=1:n1
    if mod(ii,2)==0
        array1p = [array1p C_close_p];
        array1s = [array1s C_close_s];
        array1i = [array1i C_close_i];
    else
        array1p = [array1p C_far_p];
        array1s = [array1s C_far_s];
        array1i = [array1i C_far_i];
    end
end
arrayp = [array1p C_far_p array1p C_far_p];
arrays = [array1s C_far_p array1s C_far_p];
arrayi = [array1i C_far_p array1i C_far_p];
H_p = zeros([s s]);
H_s = zeros([s s]);
H_i = zeros([s s]);

for ii=2:s
    H_p(ii, ii-1) = arrayp(ii-1);
    H_p(ii-1, ii) = arrayp(ii-1);
    H_s(ii, ii-1) = arrays(ii-1);
    H_s(ii-1, ii) = arrays(ii-1);
    H_i(ii, ii-1) = arrayi(ii-1);
    H_i(ii-1, ii) = arrayi(ii-1);
end
H_p(1,s) = C_close_p; % Connects the first and last waveguide, putting a coupling in H top right and bottom left
H_p(s,1) = C_close_p;
H_s(1,s) = C_close_s; % Connects the first and last waveguide, putting a coupling in H top right and bottom left
H_s(s,1) = C_close_s;
H_i(1,s) = C_close_i; % Connects the first and last waveguide, putting a coupling in H top right and bottom left
H_i(s,1) = C_close_i;

%*******************INTRODUCING MISSING WAVEGUIDE**************************
if missing_waveguides
    H_p(defect_position,defect_position-1) = 0; % Specify the coupling constant between the two waveguides adjacent to the defect
    H_p(defect_position-1,defect_position) = 0; % Specify the coupling constant between the two waveguides adjacent
    H_p(defect_position+1,defect_position) = 0; % Specify the coupling constant between the two waveguides adjacent to the defect
    H_p(defect_position,defect_position+1) = 0; % Specify the coupling constant between the two waveguides adjacent to the defect
    H_p(defect_position+1,defect_position-1) = 3368.598; % Specify the coupling constant between the two waveguides adjacent to the defect
    H_p(defect_position-1,defect_position+1) = 3368.598; % Specify the coupling constant between the two waveguides adjacent to the defect
    H_s(defect_position,defect_position-1) = 0; % Specify the coupling constant between the two waveguides adjacent to the defect
    H_s(defect_position-1,defect_position) = 0; % Specify the coupling constant between the two waveguides adjacent
    H_s(defect_position+1,defect_position) = 0; % Specify the coupling constant between the two waveguides adjacent to the defect
    H_s(defect_position,defect_position+1) = 0; % Specify the coupling constant between the two waveguides adjacent to the defect
    H_s(defect_position+1,defect_position-1) = 3220.895; % Specify the coupling constant between the two waveguides adjacent to the defect
    H_s(defect_position-1,defect_position+1) = 3220.895; % Specify the coupling constant between the two waveguides adjacent to the defect
    H_i(defect_position,defect_position-1) = 0; % Specify the coupling constant between the two waveguides adjacent to the defect
    H_i(defect_position-1,defect_position) = 0; % Specify the coupling constant between the two waveguides adjacent
    H_i(defect_position+1,defect_position) = 0; % Specify the coupling constant between the two waveguides adjacent to the defect
    H_i(defect_position,defect_position+1) = 0; % Specify the coupling constant between the two waveguides adjacent to the defect
    H_i(defect_position+1,defect_position-1) = 3521.412; % Specify the coupling constant between the two waveguides adjacent to the defect
    H_i(defect_position-1,defect_position+1) = 3521.412;
end

if wider_center_waveguide
    H_p((s+1)/2,(s+1)/2) =  H_p((s+1)/2,[(s+1)/2]+1)*w_factor;
    H_s((s+1)/2,(s+1)/2) =  H_s((s+1)/2,[(s+1)/2]+1)*w_factor;
    H_i((s+1)/2,(s+1)/2) =  H_i((s+1)/2,[(s+1)/2]+1)*w_factor;
end

if disorder
    for ii=2:s
        if 0
            rando = 2*rand(1)-1;
            H_p(ii, ii-1) = arrayp(ii-1)+rando*noise_lvl;
            H_p(ii-1, ii) = arrayp(ii-1)+rando*noise_lvl;
            H_s(ii, ii-1) = arrays(ii-1)+rando*noise_lvl;
            H_s(ii-1, ii) = arrays(ii-1)+rando*noise_lvl;
            H_i(ii, ii-1) = arrayi(ii-1)+rando*noise_lvl;
            H_i(ii-1, ii) = arrayi(ii-1)+rando*noise_lvl;
        end
        rando = 2*rand(1)-1;
        rando_percentage(ii)= 1+rando*percent;
        H_p(ii, ii-1) = arrayp(ii-1)*rando_percentage(ii);
        H_p(ii-1, ii) = arrayp(ii-1)*rando_percentage(ii);
        
        H_s(ii, ii-1) = arrays(ii-1)*rando_percentage(ii);
        H_s(ii-1, ii) = arrays(ii-1)*rando_percentage(ii);
        
        H_i(ii, ii-1) = arrayi(ii-1)*rando_percentage(ii);
        H_i(ii-1, ii) = arrayi(ii-1)*rando_percentage(ii);
    end
end

if on_diagonal
    for jj =1:s
        rando = 2*rand(1)-1;
        %rando_percentage(jj)= 1+rando*percent;
        H_p(jj,jj) =  H_p(jj,jj)+cc_noise_lvl*rando;
        H_s(jj,jj) =  H_s(jj,jj)+cc_noise_lvl*rando;
        H_i(jj,jj) =   H_i(jj,jj)+cc_noise_lvl*rando;
    end
end
%**************************************************************************
%*************************************************************************

%*************** Waveguides locations ************************************
cc = diag(H_p,1);
C_p=pi./(2*cc);
D_p=interp1(coupling.coupling_length, coupling.gap, C_p, 'spline');
positions       =zeros(N,1);
acc_sum         =cumsum(D_p);
positions(2:end) = acc_sum';
positions(1)=0;
if missing_waveguides
    D_p(defect_position-1)  = D_p(defect_position-3);
    D_p(defect_position)    = 0;
    D_p(defect_position+1)  = D_p(defect_position-1)+D_p(defect_position-3);
    
    D_p(defect_position2-1)  = D_p(defect_position2-3);
    D_p(defect_position2)    = 0;
    D_p(defect_position2+1)  = D_p(defect_position2-1)+D_p(defect_position2-3);
end
%*************************************************************************

% ************************** MODES & ENERGIES *****************************
% DIagonalize the Hamiltonian giving all of the eigenvectors in V and the diagonal matrix with eigenvalues on the diagonal in D
[V D]           = eig(H_p);
central         = ceil (N/2);
[val loc]       = max(abs(V(central,:)));
energies        = diag(D);
if 0
figure(999);hold on;
bar((V(1:end,end)),'b');
xlabel('x (\mum)','fontsize',18); ylabel('Amplitude (a.u.)','fontsize',14);
%xlim([positions((s+1)/2)-2e-6 positions((s+1)/2)+2e-6]);


    figure(1+15);
    subplot(3,1,1)
    scatter(positions,zeros(N,1),40,'filled');xlim([positions((s+1)/2)-2e-6 positions((s+1)/2)+2e-6]);
    relative_energy =  abs(V(51,51).^2)/(sum(abs(V(:,51)).^2)-abs(V(51,51).^2));
    subplot(3,1,2)
    energies = diag(D);
    plot(energies(1:end),'ok'); hold on;
    plot(1, energies(1),'xg');plot(s,energies(end),'xr');plot((s+1)/2,energies((s+1)/2),'xb');xlim([0 s+1]);
    xlabel('x (\mum)','fontsize',18); ylabel('Energy','fontsize',14);
    subplot(3,1,3)
    bar(positions,(V(1:end,1:end)));
    xlabel('x (\mum)','fontsize',18); ylabel('Amplitude (a.u.)','fontsize',14);xlim([positions((s+1)/2)-2e-6 positions((s+1)/2)+2e-6]);
end
%*************************************************************************

% ************************** INPUT FIELD *********************************
%input pump field
A=zeros(N,nz);
A((N+1)/2,1)=sqrt(P_in);

%input signal/idler field
B=zeros(N,N,nz);
C=zeros(N,nz);
%B((N+1)/2,(N+1)/2,1)=1; input pairs in central wg
% ************************************************************************


% ************************************************************************
% ************************** MAIN LOOP ***********************************
% ************************************************************************
for j=1:nz-1
    %B(:,:,j+1)=B(:,:,j)+1i*H_s*dz*B(:,:,j)+1i*dz*B(:,:,j)*H_i+2i*g*dz*abs(A(:,j)).^2*transpose(abs(A(:,j)).^2).*B(:,:,j)+1i*g*dz*diag(A(:,j).^2);
    %A(:,j+1)=A(:,j)+1i*H_p*dz*A(:,j)+1i*g*dz*abs(A(:,j)).^2.*A(:,j);
    A(:,j+1)=A(:,j)+1i*H_p*dz*A(:,j)-0.5*dz^2*H_p^2*A(:,j);
    B(:,:,j+1)=B(:,:,j)+1i*H_s*dz*B(:,:,j)+1i*dz*B(:,:,j)*H_i-0.5*dz^2*H_s^2*B(:,:,j)-0.5*dz^2*B(:,:,j)*H_i^2-dz^2*H_s*B(:,:,j)*H_i+1i*g*dz*diag(A(:,j).^2);
    C(:,j+1)=diag(B(:,:,j+1));
end
% ************************************************************************
% ************************************************************************
% ************************************************************************

if 0
    figure(1);
    aa=surf(0:dz:L-dz,[-(N-1)/2:(N-1)/2],abs(A).^2);xlim([0 500]*1e-6); ylim([-(s-1)/2 +(s-1)/2]);
    set(aa, 'edgecolor','none')
    view(0, 90);
    
    figure(2);A
    bb=surf([-(N-1)/2:(N-1)/2],[-(N-1)/2:(N-1)/2],abs(B(:,:,nz)).^2);%xlim([-(s-1)/2 +(s-1)/2]);ylim([-(s-1)/2 +(s-1)/2]);
    xlim([-15 15])  ; ylim([-15 15])
    set(bb, 'edgecolor','none')
    view(0, 90);
    
    figure(3);
    bb=surf(0:dz:L-dz,[-(N-1)/2:(N-1)/2],abs(C).^2);xlim([0 500]*1e-6); %ylim([-(s-1)/2 +(s-1)/2]);
    ylim([-15 15])
    set(bb, 'edgecolor','none')
    view(0, 90);
end

figure(1000);
fig = gcf;
set(fig,'PaperUnits','inches');
set(fig,'PaperPosition',[0 0 10 6]);
subplot(2,3,1);
scatter(positions,zeros(N,1),40,'filled');xlim([positions((s+1)/2)-1e-6 positions((s+1)/2)+1e-6]); xlabel('x (m)');title('Waveguides locations');
subplot(2,3,4);
plot(energies(1:end),'ok'); hold on;plot(1, energies(1),'xg');plot(s,energies(end),'xr');plot((s+1)/2,energies((s+1)/2),'xb');xlim([0 s+1]);xlabel('Mode number'); ylabel('k_x (m^-^1)');title('Modes energies');
subplot(2,3,[2 3]);aa=surf(0:dz:L-dz,[-(N-1)/2:(N-1)/2],abs(A).^2);xlim([0 500]*1e-6); ylim([-(s-1)/2 +(s-1)/2]); ylim([-(s-1)/2 +(s-1)/2]);set(aa, 'edgecolor','none');view(0, 90);xlabel('z (m)'); ylabel('Waveguide number');title('Pump power propagation');
subplot(2,3,[5 6]);bb=surf(0:dz:L-dz,[-(N-1)/2:(N-1)/2],abs(C).^2);xlim([0 500]*1e-6); ylim([-(s-1)/2 +(s-1)/2]);set(bb, 'edgecolor','none');view(0, 90);xlabel('z (m)'); ylabel('Waveguide number');title('Photon pairs propagation');
if SAVEFIGS
    saveas(fig,'BGnProp_2_4_trivial_ondiagonal1_5', 'tiff');
    saveas(fig,'BGnProp_2_4_trivial_ondiagonal1_5', 'fig');
    %saveas(fig,'2_7_topo_noise0e4', 'tiff');
    %saveas(fig,'2_7_topo_noise0e4', 'fig');
end
figure(1001);
bb=surf([-(N-1)/2:(N-1)/2],[-(N-1)/2:(N-1)/2],abs(B(:,:,nz)).^2);%xlim([-(s-1)/2 +(s-1)/2]);ylim([-(s-1)/2 +(s-1)/2]);
xlim([-10 10])  ; ylim([-10 10]); 
caxis([0 4.8e-4]);
set(bb, 'edgecolor','none')
view(0, 90);
if SAVEFIGS
    saveas(bb,'Co_2_4_trivial_ondiagonal1_5', 'tiff');
    saveas(bb,'Co_2_4_trivial_ondiagonal1_5', 'fig');
    %saveas(bb,'Co_2_7_topo_noise0e4', 'tiff');
    %saveas(bb,'Co_2_7_topo_noise0e4', 'fig');
end

figure(1002);% Plots the intensity distribution at the output
%fig1002 = plot((abs(C(98:106, end)).^2)/[max(abs(C(98:106, end)).^2)]);
fig1002 = plot(abs(C(98:106, end)).^2);
if SAVEFIGS
    saveas(fig1002,'Out_2_7_topo_ondiagonal2', 'tiff');
    saveas(fig1002,'Out_2_7_topo_ondiagonal2', 'fig');
    %saveas(fig1002,'Co_2_7_topo_noise0e4', 'tiff');
    %saveas(fig1002,'Co_2_7_topo_noise0e4', 'fig');
end


