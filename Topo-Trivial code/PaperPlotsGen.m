%close all;

rng('shuffle');

per=0;
percent=per/100;

% N is odd
N=203;% the size of lattice 
halfchain=floor(N-1)/2;

% Coupling Parameters for Hamiltonian;
%pump
    tp_s=46118;
    tp_l=14951;
%signal 
    ts_s=42882;  %%%%% LINE 291
    ts_l=13603;  %%%%% LINE 292
%idle
    ti_s=45162;  %%%%% LINE 351
    ti_l=14562;  %%%%% LINE 352
    
    g=120;


    
    
    
% Define Hamiltonian;
t1=zeros(1,halfchain-1);
t2=zeros(1,halfchain);
% generat a random number in [-1,1];   1+(2*rand(1)-1)*percent;
disorder1=1+(2*rand(1,length(t1))-1)*percent;
disorder2=1+(2*rand(1,length(t2))-1)*percent;

% disordered componets for pump


 
t1(1:2:end)=tp_s*disorder1(1:2:length(t1));
t1(2:2:end)=tp_l*disorder1(2:2:length(t1));


t2(1:2:end)=tp_s*disorder2(1:2:length(t2));
t2(2:2:end)=tp_l*disorder2(2:2:length(t2));
%t2(1:2:end)=tp_l*(1+(2*rand(1,length(1:2:length(t2)))-1)*percent);
%t2(2:2:end)=tp_s*(1+(2*rand(1,length(2:2:length(t2)))-1)*percent);

dis_mid=1+(2*rand-1)*percent;
tp=[t1,tp_s*dis_mid,t2];

dis_end=1+(2*rand-1)*percent;
Hp=diag(tp,1);
Hp(end:1)=tp_l*dis_end; % periodic boundary
Hp=Hp+Hp';
[Vp,Dp]=eig(Hp);


% componets for signal
t1(1:2:end)=ts_s*disorder1(1:2:length(t1));
t1(2:2:end)=ts_l*disorder1(2:2:length(t1));

t2(1:2:end)=ts_s*disorder2(1:2:length(t2));
t2(2:2:end)=ts_l*disorder2(2:2:length(t2));

ts=[t1,ts_s*dis_mid,t2];

Hs=diag(ts,1);
Hs(end:1)=ts_l*dis_end; % periodic boundary
Hs=Hs+Hs';
[Vs,Ds]=eig(Hs);

% componets for idler
t1(1:2:end)=ti_s*disorder1(1:2:length(t1));
t1(2:2:end)=ti_l*disorder1(2:2:length(t1));

t2(1:2:end)=ti_s*disorder2(1:2:length(t2));
t2(2:2:end)=ti_l*disorder2(2:2:length(t2));

ti=[t1,ti_s*dis_mid,t2];

% short-short defect;
Hi= diag(ti,1);
Hi(end:1)=ti_l*dis_end; % periodic boundary
Hi=Hi+Hi';
[Vi,Di]=eig(Hi);

Tp_p=Vp(:,halfchain+1);

Tp_s=Vs(:,halfchain+1);
Tr1_s=Vs(:,1);
Tr2_s=Vs(:,end);

Tp_i=Vi(:,halfchain+1);
Tr1_i=Vi(:,1);
Tr2_i=Vi(:,end);

 TpTp_si=Tp_s*Tp_i.'; %f1
 TpTr1_si=Tp_s*Tr1_i.'; %f6
 TpTr2_si=Tp_s*Tr2_i.'; %f7 
 Tr1Tp_si=Tr1_s*Tp_i.'; %f8
 Tr2Tp_si=Tr2_s*Tp_i.'; %f9
 
 Tr1Tr1_si=Tr1_s*Tr1_i.'; %f2
 Tr1Tr2_si=Tr1_s*Tr2_i.'; %f3
 Tr2Tr1_si=Tr2_s*Tr1_i.'; %f4
 Tr2Tr2_si=Tr2_s*Tr2_i.'; %f5
 
 TrTr_si(:,:,1)=TpTp_si;
 TrTr_si(:,:,2)=TpTr1_si;
 TrTr_si(:,:,3)=TpTr2_si;
 TrTr_si(:,:,4)=Tr1Tp_si;
 TrTr_si(:,:,5)=Tr2Tp_si;
 
 TrTr_si(:,:,6)=Tr1Tr1_si;
 TrTr_si(:,:,7)=Tr1Tr2_si;
 TrTr_si(:,:,8)=Tr2Tr1_si;
 TrTr_si(:,:,9)=Tr2Tr2_si; 

%TrTr_si=zeros(N,N,N^2);

% for i = 1:N
%     for j=1:N
% 
%     TrTr_si(:,:,(i-1)*N+j)=Vs(:,i)*Vi(:,j).';
% 
%     end
% end
% the evolution dynamics 
% parametres from Andrea;
L=0.381E-3; %length of structure
%L=0.003E-3;
nz = 1E3; %number of steps along length
dt=L/nz; %size of a step

time=dt*(1:nz);




   
% define the input sates

psi_p=zeros(N,length(time));
psi_p(halfchain+1,1)=1;% halfchain is the adjacent position; halfchain+1 is the center of the short-short defect

psi_si=zeros(N,N,length(time));

Us=Vs*diag(exp(-1i*diag(Ds)*dt))*Vs^(-1);
Ui=Vi*diag(exp(-1i*diag(Di)*dt))*Vi^(-1);
Up=Vp*diag(exp(-1i*diag(Dp)*dt))*Vp^(-1);
%figure;imagesc([85:120],[85:120],abs(psi_si(85:120,85:120,1)));mm(:,1)=getframe(gcf);

ff=zeros(9,length(time));
SS=zeros(length(time),1);

for T=2:length(time)
    
 psi_p(:,T)=Up*psi_p(:,T-1);
% H=Hs+Hi+g*diag(psi_p(:,time).^2);

  psi_si(:,:,T)=Us*psi_si(:,:,T-1)*Ui+diag(1i*g*dt*psi_p(:,T-1).^2);
% psi_si(:,:,T)=(Us-eye(N))*psi_si(:,:,T-1)+psi_si(:,:,T-1)*(Ui)-dt^2*Hs*psi_si(:,:,T-1)*Hi+1i*g*dt*diag(psi_p(:,T-1).^2);
% pc(:,T)=diag(psi_si(:,:,T));
% imagesc([100:104],[100:104],abs(psi_si(100:104,100:104,T)));title(['time is ',num2str(T)]);colorbar;pause(0);mm(:,T)=getframe(gcf);
%pause(0.1);
  for k=1:9
  
      ff(k,T)=sum(sum(conj(psi_si(:,:,T)) .* TrTr_si(:,:,k)));
      
  end
   SS(T)=sum(sum(abs(psi_si(:,:,T)).^2));
end



%%
psi_p60=psi_p;
psi_si60=psi_si;
Tp_p60=Tp_p;
Tr1_p60=Vp(:,1);
Tr2_p60=Vp(:,end);
ff60=ff;
Dp60=Dp;
%%
psi_p40=psi_p;
psi_si40=psi_si;
Tp_p40=Tp_p;
Tr1_p40=Vp(:,1);
Tr2_p40=Vp(:,end);
ff40=ff;
Dp40=Dp;
%%
psi_p20=psi_p;
psi_si20=psi_si;
Tp_p20=Tp_p;
Tr1_p20=Vp(:,1);
Tr2_p20=Vp(:,end);
ff20=ff;
Dp20=Dp;
%% Fig.1 

figure('position', [0, 0, 1600, 800]);title(['Dyanmics of Pump and correlated biphotos'])
subplot(2,1,1);imagesc(abs(psi_p(94:110,1:1000)).^2);colormap(hot);colorbar;
title(['Dyanmics of Pump ',num2str(per),' percent disorder']);
yticks([7 8 9 10 11])
yticklabels({'100','101','102','103','104'})
xticks([5 920 975]);
xticklabels({'A','B','C'});

hold on;
X=5*ones(1,203);
Y=1:203;
plot(X,Y,'w-.','LineWidth',3);set(gca,'fontsize', 13);

X=920*ones(1,203);
Y=1:203;
plot(X,Y,'w-.','LineWidth',3);set(gca,'fontsize', 13);

X=975*ones(1,203);
Y=1:203;
plot(X,Y,'w-.','LineWidth',3);set(gca,'fontsize', 13);

for time=1:1000
    
psi_0(:,time)=sum(abs(psi_si(:,:,time)).^2);
psi_1(:,time)=diag(abs(psi_si(:,:,time))).^2;


end

%subplot(3,1,2);imagesc(abs(psi_0(94:110,1:1000)));colormap(hot);colorbar;
subplot(2,1,2);imagesc(abs(psi_1(94:110,1:1000)));colormap(hot);colorbar;
title(['Dyanmics of correlated biphotos'])
yticks([7 8 9 10 11])
yticklabels({'100','101','102','103','104'})

hold on;
X=5*ones(1,203);
Y=1:203;
plot(X,Y,'w-.','LineWidth',3);set(gca,'fontsize', 13);
X=920*ones(1,203);
Y=1:203;
plot(X,Y,'w-.','LineWidth',3);set(gca,'fontsize', 13);
X=975*ones(1,203);
Y=1:203;
plot(X,Y,'w-.','LineWidth',3);set(gca,'fontsize', 13);

xticks([5 920 975]);
xticklabels({'A','B','C'});set(gca,'fontsize', 13);
%% Fig.2
%close all;
figure('position', [0, 0, 1600, 800]);title(['tensor of Pump at different step'])
l=101:103;
T=5;
P_5=kron(psi_p(:,T),psi_p(:,T)');
subplot(2,3,1);imagesc(abs(P_5(l,l)).^2);colormap(jet);colorbar;title(['A of pump']);set(gca,'fontsize', 13);
xticks([1 2 3]);
xticklabels({'-1','0','1'});
yticks([1 2 3]);
yticklabels({'-1','0','1'});
xlabel(['n_p']);ylabel(['n_p'])
T=920;
P_5=kron(psi_p(:,T),psi_p(:,T)');
subplot(2,3,2);imagesc(abs(P_5(l,l)).^2);colormap(jet);colorbar;title(['B of pump']);set(gca,'fontsize', 13);
xticks([1 2 3]);
xticklabels({'-1','0','1'});
yticks([1 2 3]);
yticklabels({'-1','0','1'});
T=975;
P_5=kron(psi_p(:,T),psi_p(:,T)');
subplot(2,3,3);imagesc(abs(P_5(l,l)).^2);colormap(jet);colorbar;title(['C of pump']);set(gca,'fontsize', 13);
xticks([1 2 3]);
xticklabels({'-1','0','1'});
yticks([1 2 3]);
yticklabels({'-1','0','1'});


T=5;
P_5=psi_si(:,:,T);
subplot(2,3,4);imagesc(abs(P_5(l,l)).^2);colormap(jet);colorbar;title(['A of biphoton']);set(gca,'fontsize', 13);
xticks([1 2 3]);
xticklabels({'-1','0','1'});
yticks([1 2 3]);
yticklabels({'-1','0','1'});
xlabel(['n_s']);ylabel(['n_i'])
T=920;
P_5=psi_si(:,:,T);
subplot(2,3,5);imagesc(abs(P_5(l,l)).^2);colormap(jet);colorbar;title(['B of biphoton']);set(gca,'fontsize', 13);
xticks([1 2 3]);
xticklabels({'-1','0','1'});
yticks([1 2 3]);
yticklabels({'-1','0','1'});
T=975;
P_5=psi_si(:,:,T);
subplot(2,3,6);imagesc(abs(P_5(l,l)).^2);colormap(jet);colorbar;title(['C of biphoton']);set(gca,'fontsize', 13);
xticks([1 2 3]);
xticklabels({'-1','0','1'});
yticks([1 2 3]);
yticklabels({'-1','0','1'});

%% Fig.3 
%figure('position', [0, 0, 1600, 800]);title(['tensor of Pump at different step'])

TT=[5,920,975];
F=zeros(3,3);
FF_si=zeros(3,3,3);
FF_bump=zeros(3,3);
Tr1_p=Vp(:,1);
Tr2_p=Vp(:,end);

for i=1 :3
  T=TT(i);
  F(1,1)= abs(ff(1,T)).^2; 
  F(1,2)= abs(ff(2,T)).^2; 
  F(1,3)= abs(ff(3,T)).^2; 
  F(2,1)= abs(ff(4,T)).^2;
  F(2,2)= abs(ff(6,T)).^2; 
  F(2,3)= abs(ff(7,T)).^2;
  F(3,1)= abs(ff(5,T)).^2;
  F(3,2)= abs(ff(8,T)).^2;
  F(3,3)= abs(ff(9,T)).^2;
  FF_si(:,:,i)=F;
  FF_bump(1,i)=abs(sum(conj(psi_p(:,T)).*Tp_p)).^2;
  FF_bump(2,i)=abs(sum(conj(psi_p(:,T)).*Tr1_p)).^2;
  FF_bump(3,i)=abs(sum(conj(psi_p(:,T)).*Tr2_p)).^2;
end



figure('position', [0, 0, 1600, 800]);title(['tensor of Pump at different step'])
subplot(2,3,1);bar(FF_bump(:,1));title(['A of Pump']);xticklabels({'Tp','Tr1','Tr2'});set(gca,'fontsize', 13);
subplot(2,3,2);bar(FF_bump(:,2));title(['B of Pump']);xticklabels({'Tp','Tr1','Tr2'});set(gca,'fontsize', 13);
subplot(2,3,3);bar(FF_bump(:,3));title(['C of Pump']);xticklabels({'Tp','Tr1','Tr2'});set(gca,'fontsize', 13);



% for biphoton states
subplot(2,3,4);imagesc(FF_si(:,:,1));colormap(jet);colorbar;title(['A of biphoton']);
xticks([1 2 3]);yticks([1 2 3]);set(gca,'fontsize', 13);
xticklabels({'Tp','Tr1','Tr2'});yticklabels({'Tp','Tr1','Tr2'})
subplot(2,3,5);imagesc(FF_si(:,:,2));colormap(jet);colorbar;title(['B of biphoton']);
xticks([1 2 3]);yticks([1 2 3]);set(gca,'fontsize', 13);
xticklabels({'Tp','Tr1','Tr2'});yticklabels({'Tp','Tr1','Tr2'})
subplot(2,3,6);imagesc(FF_si(:,:,3));colormap(jet);colorbar;title(['C of biphoton']);
xticks([1 2 3]);yticks([1 2 3]);set(gca,'fontsize', 13);
xticklabels({'Tp','Tr1','Tr2'});yticklabels({'Tp','Tr1','Tr2'})

%% Fig 4. disorder 20%, 40%, 60% WaveGuideModes;  29 OCT 2019

T=975;
psi_p=psi_p20;
psi_si=psi_si20;

figure('position', [0, 0, 1600, 800]);
P_5=kron(psi_p(:,T),psi_p(:,T)');
subplot(2,3,1);imagesc(abs(P_5(l,l)).^2);colormap(jet);colorbar;title(['C of pump--20% disorder']);set(gca,'fontsize', 13);
xticks([1 2 3]);
xticklabels({'-1','0','1'});
yticks([1 2 3]);
yticklabels({'-1','0','1'});
xlabel(['n_p']);ylabel(['n_p'])

P_5=psi_si(:,:,T);
subplot(2,3,4);imagesc(abs(P_5(l,l)).^2);colormap(jet);colorbar;title(['C of biphoton--20% disorder']);set(gca,'fontsize', 13);
xticks([1 2 3]);
xticklabels({'-1','0','1'});
yticks([1 2 3]);
yticklabels({'-1','0','1'});
xlabel(['n_s']);ylabel(['n_i'])


psi_p=psi_p40;
psi_si=psi_si40;


P_5=kron(psi_p(:,T),psi_p(:,T)');
subplot(2,3,2);imagesc(abs(P_5(l,l)).^2);colormap(jet);colorbar;title(['C of pump--40% disorder']);set(gca,'fontsize', 13);
xticks([1 2 3]);
xticklabels({'-1','0','1'});
yticks([1 2 3]);
yticklabels({'-1','0','1'});
xlabel(['n_p']);ylabel(['n_p'])

P_5=psi_si(:,:,T);
subplot(2,3,5);imagesc(abs(P_5(l,l)).^2);colormap(jet);colorbar;title(['C of biphoton--40% disorder']);set(gca,'fontsize', 13);
xticks([1 2 3]);
xticklabels({'-1','0','1'});
yticks([1 2 3]);
yticklabels({'-1','0','1'});
xlabel(['n_s']);ylabel(['n_i'])

psi_p=psi_p60;
psi_si=psi_si60;


P_5=kron(psi_p(:,T),psi_p(:,T)');
subplot(2,3,3);imagesc(abs(P_5(l,l)).^2);colormap(jet);colorbar;title(['C of pump--60% disorder']);set(gca,'fontsize', 13);
xticks([1 2 3]);
xticklabels({'-1','0','1'});
yticks([1 2 3]);
yticklabels({'-1','0','1'});
xlabel(['n_p']);ylabel(['n_p'])

P_5=psi_si(:,:,T);
subplot(2,3,6);imagesc(abs(P_5(l,l)).^2);colormap(jet);colorbar;title(['C of biphoton--60% disorder']);set(gca,'fontsize', 13);
xticks([1 2 3]);
xticklabels({'-1','0','1'});
yticks([1 2 3]);
yticklabels({'-1','0','1'});
xlabel(['n_s']);ylabel(['n_i'])

%%  Fig 5. disorder 20%, 40%, 60% EigenModes;  29 OCT 2019

T=975;

psi_p=psi_p20;
psi_si=psi_si20;
Tp_p=Tp_p20;
Tr1_p=Tr1_p20;
Tr2_p=Tr2_p20;
ff=ff20;

FF_si=zeros(3,3);
FF_bump=zeros(3);

  FF_si(1,1)= abs(ff(1,T)).^2; 
  FF_si(1,2)= abs(ff(2,T)).^2; 
  FF_si(1,3)= abs(ff(3,T)).^2; 
  FF_si(2,1)= abs(ff(4,T)).^2;
  FF_si(2,2)= abs(ff(6,T)).^2; 
  FF_si(2,3)= abs(ff(7,T)).^2;
  FF_si(3,1)= abs(ff(5,T)).^2;
  FF_si(3,2)= abs(ff(8,T)).^2;
  FF_si(3,3)= abs(ff(9,T)).^2;

  FF_bump(1)=abs(sum(conj(psi_p(:,T)).*Tp_p)).^2;
  FF_bump(2)=abs(sum(conj(psi_p(:,T)).*Tr1_p)).^2;
  FF_bump(3)=abs(sum(conj(psi_p(:,T)).*Tr2_p)).^2;

figure('position', [0, 0, 1600, 800]);title(['tensor of Pump at different step'])

subplot(2,3,1);bar(FF_bump);title(['C of Pump--disorder 20%']);xticklabels({'Tp','Tr1','Tr2'});set(gca,'fontsize', 13);

% for biphoton states
subplot(2,3,4);imagesc(FF_si(:,:));colormap(jet);colorbar;title(['C of biphoton--disorder 20%']);
xticks([1 2 3]);yticks([1 2 3]);set(gca,'fontsize', 13);
xticklabels({'Tp','Tr1','Tr2'});yticklabels({'Tp','Tr1','Tr2'})



psi_p=psi_p40;
psi_si=psi_si40;
Tp_p=Tp_p40;
Tr1_p=Tr1_p40;
Tr2_p=Tr2_p40;
ff=ff40;

FF_si=zeros(3,3);
FF_bump=zeros(3);

  FF_si(1,1)= abs(ff(1,T)).^2; 
  FF_si(1,2)= abs(ff(2,T)).^2; 
  FF_si(1,3)= abs(ff(3,T)).^2; 
  FF_si(2,1)= abs(ff(4,T)).^2;
  FF_si(2,2)= abs(ff(6,T)).^2; 
  FF_si(2,3)= abs(ff(7,T)).^2;
  FF_si(3,1)= abs(ff(5,T)).^2;
  FF_si(3,2)= abs(ff(8,T)).^2;
  FF_si(3,3)= abs(ff(9,T)).^2;

  FF_bump(1)=abs(sum(conj(psi_p(:,T)).*Tp_p)).^2;
  FF_bump(2)=abs(sum(conj(psi_p(:,T)).*Tr1_p)).^2;
  FF_bump(3)=abs(sum(conj(psi_p(:,T)).*Tr2_p)).^2;


subplot(2,3,2);bar(FF_bump);title(['C of Pump--disorder 40%']);xticklabels({'Tp','Tr1','Tr2'});set(gca,'fontsize', 13);

% for biphoton states
subplot(2,3,5);imagesc(FF_si(:,:));colormap(jet);colorbar;title(['C of biphoton--disorder 40%']);
xticks([1 2 3]);yticks([1 2 3]);set(gca,'fontsize', 13);
xticklabels({'Tp','Tr1','Tr2'});yticklabels({'Tp','Tr1','Tr2'})



psi_p=psi_p60;
psi_si=psi_si60;
Tp_p=Tp_p60;
Tr1_p=Tr1_p60;
Tr2_p=Tr2_p60;
ff=ff60;

FF_si=zeros(3,3);
FF_bump=zeros(3);

  FF_si(1,1)= abs(ff(1,T)).^2; 
  FF_si(1,2)= abs(ff(2,T)).^2; 
  FF_si(1,3)= abs(ff(3,T)).^2; 
  FF_si(2,1)= abs(ff(4,T)).^2;
  FF_si(2,2)= abs(ff(6,T)).^2; 
  FF_si(2,3)= abs(ff(7,T)).^2;
  FF_si(3,1)= abs(ff(5,T)).^2;
  FF_si(3,2)= abs(ff(8,T)).^2;
  FF_si(3,3)= abs(ff(9,T)).^2;

  FF_bump(1)=abs(sum(conj(psi_p(:,T)).*Tp_p)).^2;
  FF_bump(2)=abs(sum(conj(psi_p(:,T)).*Tr1_p)).^2;
  FF_bump(3)=abs(sum(conj(psi_p(:,T)).*Tr2_p)).^2;


subplot(2,3,3);bar(FF_bump);title(['C of Pump--disorder 60%']);xticklabels({'Tp','Tr1','Tr2'});set(gca,'fontsize', 13);

% for biphoton states
subplot(2,3,6);imagesc(FF_si(:,:));colormap(jet);colorbar;title(['C of biphoton--disorder 60%']);
xticks([1 2 3]);yticks([1 2 3]);set(gca,'fontsize', 13);
xticklabels({'Tp','Tr1','Tr2'});yticklabels({'Tp','Tr1','Tr2'})

%% Fig.6 BandStructure and Disorders
% 24 Sep 2019; the plot of band structure while increasing the disorder
% strength. and using the red color to mark the eigenmodes contributing the
% generation of TPTP biphoton states


N=203;% the size of lattice 
halfchain=floor(N-1)/2;

percent0=[0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.6,0.7,0.8,0.9,1];

% Coupling Parameters for Hamiltonian;
%pump
    tp_s=46118;
    tp_l=14951;
%signal 
    ts_s=42882;  %%%%% LINE 291
    ts_l=13603;  %%%%% LINE 292
%idle
    ti_s=45162;  %%%%% LINE 351
    ti_l=14562;  %%%%% LINE 352
    
    g=120;
DDp=zeros(N,length(percent0));
DDs=zeros(N,length(percent0));
DDi=zeros(N,length(percent0));
Ptr_p=zeros(N,length(percent0));
Ptr_s=zeros(N,length(percent0));
Ptr_i=zeros(N,length(percent0));

for i =1:length(percent0)    
 
 percent=percent0(i);
% Define Hamiltonian;
t1=zeros(1,halfchain-1);
t2=zeros(1,halfchain);
% generat a random number in [-1,1];   1+(2*rand(1)-1)*percent;
disorder1=1+(2*rand(1,length(t1))-1)*percent;
disorder2=1+(2*rand(1,length(t2))-1)*percent;

% disordered componets for pump


 
t1(1:2:end)=tp_s*disorder1(1:2:length(t1));
t1(2:2:end)=tp_l*disorder1(2:2:length(t1));


t2(1:2:end)=tp_s*disorder2(1:2:length(t2));
t2(2:2:end)=tp_l*disorder2(2:2:length(t2));
%t2(1:2:end)=tp_l*(1+(2*rand(1,length(1:2:length(t2)))-1)*percent);
%t2(2:2:end)=tp_s*(1+(2*rand(1,length(2:2:length(t2)))-1)*percent);

dis_mid=1+(2*rand-1)*percent;
tp=[t1,tp_s*dis_mid,t2];

dis_end=1+(2*rand-1)*percent;
Hp=diag(tp,1);
Hp(end:1)=tp_l*dis_end; % periodic boundary
Hp=Hp+Hp';
[Vp,Dp]=eig(Hp);

DDp(:,i)=diag(Dp);


 tr_index1=1:101;  
 tr_index2=103:203; 
 Ptr_p(:,i)=[sum((conj(Vp(:,tr_index1).^2).*kron(ones(1,length(tr_index1)),Vp(:,102))).^2./kron(ones(203,1),sum(abs(conj(Vp(:,tr_index1).^2)).^2))),0,...
             sum((conj(Vp(:,tr_index2).^2).*kron(ones(1,length(tr_index2)),Vp(:,102))).^2./kron(ones(203,1),sum(abs(conj(Vp(:,tr_index2).^2)).^2)))];





% componets for signal
t1(1:2:end)=ts_s*disorder1(1:2:length(t1));
t1(2:2:end)=ts_l*disorder1(2:2:length(t1));

t2(1:2:end)=ts_s*disorder2(1:2:length(t2));
t2(2:2:end)=ts_l*disorder2(2:2:length(t2));

ts=[t1,ts_s*dis_mid,t2];

Hs=diag(ts,1);
Hs(end:1)=ts_l*dis_end; % periodic boundary
Hs=Hs+Hs';
[Vs,Ds]=eig(Hs);

 DDs(:,i)=diag(Ds);
 
 
 Ptr_s(:,i)=[sum((conj(Vs(:,tr_index1).^2).*kron(ones(1,length(tr_index1)),Vs(:,102))).^2./kron(ones(203,1),sum(abs(conj(Vs(:,tr_index1).^2)).^2))),0,...
             sum((conj(Vs(:,tr_index2).^2).*kron(ones(1,length(tr_index2)),Vs(:,102))).^2./kron(ones(203,1),sum(abs(conj(Vs(:,tr_index2).^2)).^2)))];


% componets for idler
t1(1:2:end)=ti_s*disorder1(1:2:length(t1));
t1(2:2:end)=ti_l*disorder1(2:2:length(t1));

t2(1:2:end)=ti_s*disorder2(1:2:length(t2));
t2(2:2:end)=ti_l*disorder2(2:2:length(t2));

ti=[t1,ti_s*dis_mid,t2];

% short-short defect;
Hi=diag(ti,1);
Hi(end:1)=ti_l*dis_end; % periodic boundary
Hi=Hi+Hi';
[Vi,Di]=eig(Hi);

DDi(:,i)=diag(Di);


Ptr_i(:,i)=[sum((conj(Vi(:,tr_index1).^2).*kron(ones(1,length(tr_index1)),Vi(:,102))).^2./kron(ones(203,1),sum(abs(conj(Vi(:,tr_index1).^2)).^2))),0,...
            sum((conj(Vi(:,tr_index2).^2).*kron(ones(1,length(tr_index2)),Vi(:,102))).^2./kron(ones(203,1),sum(abs(conj(Vi(:,tr_index2).^2)).^2)))];

end
close all;

%figure('position', [50, 60, 500, 300]);


for i = 1: length(percent0)
plot(percent0(i)*ones(203,1),DDp(:,i),'.'); hold on;
for j =1: N
    if Ptr_p(j,i)>0.025
     plot(percent0(i),DDp(j,i),'p','Color','r','MarkerSize',8);  
    end
end

end
xlabel('disorder');ylabel('Energy of pump');

set(gca,'fontsize', 13);

%% Fig.7 Probability of each Eigen components; 29 OCT 2019;
close all;

figure('position', [0, 0, 1200, 1000]);
l=100:104;

subplot(3,3,1);set(gca,'fontsize', 13);imagesc(abs(TpTp_si(l,l)).^2);colormap(jet);colorbar;ylabel('n_i');
xticks([1 2 3 4 5]);
xticklabels({'-2','-1','0','1', '2'});
yticks([1 2 3 4 5]);
yticklabels({'-2','-1','0','1', '2'});
subplot(3,3,2);set(gca,'fontsize', 13);imagesc(abs(TpTr1_si(l,l)).^2);colormap(jet);colorbar;
xticks([1 2 3 4 5]);
xticklabels({'-2','-1','0','1', '2'});
yticks([1 2 3 4 5]);
yticklabels({'-2','-1','0','1', '2'});
subplot(3,3,3);set(gca,'fontsize', 13);imagesc(abs(TpTr2_si(l,l)).^2);colormap(jet);colorbar;
xticks([1 2 3 4 5]);
xticklabels({'-2','-1','0','1', '2'});
yticks([1 2 3 4 5]);
yticklabels({'-2','-1','0','1', '2'});
subplot(3,3,4);set(gca,'fontsize', 13);imagesc(abs(Tr1Tp_si(l,l)).^2);colormap(jet);colorbar;ylabel('n_i');
xticks([1 2 3 4 5]);
xticklabels({'-2','-1','0','1', '2'});
yticks([1 2 3 4 5]);
yticklabels({'-2','-1','0','1', '2'});
subplot(3,3,5);set(gca,'fontsize', 13);imagesc(abs(Tr1Tr1_si(l,l)).^2);colormap(jet);colorbar;
xticks([1 2 3 4 5]);
xticklabels({'-2','-1','0','1', '2'});
yticks([1 2 3 4 5]);
yticklabels({'-2','-1','0','1', '2'});
subplot(3,3,6);set(gca,'fontsize', 13);imagesc(abs(Tr1Tr2_si(l,l)).^2);colormap(jet);colorbar;
xticks([1 2 3 4 5]);
xticklabels({'-2','-1','0','1', '2'});
yticks([1 2 3 4 5]);
yticklabels({'-2','-1','0','1', '2'});
subplot(3,3,7);set(gca,'fontsize', 13);imagesc(abs(Tr2Tp_si(l,l)).^2);colormap(jet);colorbar;ylabel('n_i');xlabel('n_s');
xticks([1 2 3 4 5]);
xticklabels({'-2','-1','0','1', '2'});
yticks([1 2 3 4 5]);
yticklabels({'-2','-1','0','1', '2'});
subplot(3,3,8);set(gca,'fontsize', 13);imagesc(abs(Tr2Tr1_si(l,l)).^2);colormap(jet);colorbar;xlabel('n_s');
xticks([1 2 3 4 5]);
xticklabels({'-2','-1','0','1', '2'});
yticks([1 2 3 4 5]);
yticklabels({'-2','-1','0','1', '2'});
subplot(3,3,9);set(gca,'fontsize', 13);imagesc(abs(Tr2Tr2_si(l,l)).^2);colormap(jet);colorbar;xlabel('n_s');
xticks([1 2 3 4 5]);
xticklabels({'-2','-1','0','1', '2'});
yticks([1 2 3 4 5]);
yticklabels({'-2','-1','0','1', '2'});


 [ax1,h1]=suplabel('Tp','x',[-.21 .08 .84 .84]);
 [ax1,h1]=suplabel('Tr1','x',[.08 .08 .84 .84]);
 [ax1,h1]=suplabel('Tr2','x',[.36 .08 .84 .84]);
 [ax2,h2]=suplabel('Tr2','y',[.10 -.20 .84 .84]);
 [ax2,h2]=suplabel('Tr1','y',[.10 .10 .84 .84]);
 [ax2,h2]=suplabel('Tp','y', [.10 .40 .84 .84]);
 
 writematrix(abs(TpTp_si(l,l)), 'tptp.csv')
 writematrix(abs(Tr1Tr1_si(l,l)), 'tr1tr1.csv')
 writematrix(abs(Tr2Tr2_si(l,l)), 'tr2tr2.csv')
 writematrix(abs(Tr1Tr2_si(l,l)), 'tr1tr2.csv')
 writematrix(abs(Tr2Tr1_si(l,l)), 'tr2tr1.csv')
 writematrix(abs(TpTr1_si(l,l)), 'tptr1.csv')
 writematrix(abs(TpTr2_si(l,l)), 'tptr2.csv')
 writematrix(abs(Tr1Tp_si(l,l)), 'tr1tp.csv')
 writematrix(abs(Tr2Tp_si(l,l)), 'tr2tp.csv')
%% Fig.8 Correlation between the weights Probabilities that the biphoton states 
%  are in the 9 Eigenmodes of linear Hamiltonian as shown in Fig.7, and disorder level

load('/Users/vivien/Dropbox (Sydney Uni)/Redondo_Zhang_Bartlett/HPC/Data/plot_all.mat')

close all;
 percent0=[0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50];
 Prob_end=[Prob_5percent(:,end),Prob_10percent(:,end),Prob_15percent(:,end),Prob_20percent(:,end),Prob_25percent(:,end)...
     Prob_30percent(:,end),Prob_35percent(:,end),Prob_40percent(:,end),Prob_45percent(:,end),Prob_50percent(:,end)];
 plot(percent0,Prob_end(1,:),'o-');hold on;
  plot(percent0,Prob_end(7,:),'.','MarkerSize',14);hold on;
  plot(percent0,Prob_end(8,:),'d-.','MarkerSize',8);hold on;

  plot(percent0,Prob_end(9,:),'>--','MarkerSize',11);hold on;
 
  plot(percent0,Prob_end(6,:),'*');hold on;
   
  plot(percent0,Prob_end(3,:),'.k--','MarkerSize',8);hold on;
  
 plot(percent0,Prob_end(5,:),'o');hold on;
 
 
 plot(percent0,Prob_end(4,:),'o','MarkerSize',10);hold on;
 plot(percent0,Prob_end(2,:),'o','MarkerSize',15);hold on;

xlim([0.03 0.53]);
ylim([-0.03 0.32]);

xlabel('disorder');ylabel('Probability of different eigenmodes');
legend('TpTp','Tr1Tr2','Tr2Tr1','Tr2Tr2','Tr1Tr1','TpTr2','Tr2Tp','Tr1Tp','TpTr1');

set(gca,'fontsize', 13);




