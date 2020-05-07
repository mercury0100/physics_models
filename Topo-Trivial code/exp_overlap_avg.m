clear all; close all; format short

n = 2;          % Number of correlations to view either side of centre (even)
N = 203;        % Total number of waveguides (odd!)
yvals = (N+1)/2-n:(N+1)/2+n;

% Load experimental data
T2 = readtable(['~/Dropbox/Directory/Quantum/DATA/Workbook2F.xlsx']);
T3 = readtable(['~/Dropbox/Directory/Quantum/DATA/Workbook3F.xlsx']);
T4 = readtable(['~/Dropbox/Directory/Quantum/DATA/Workbook4F.xlsx']);
T5 = readtable(['~/Dropbox/Directory/Quantum/DATA/Workbook5F.xlsx']);
T7 = readtable(['~/Dropbox/Directory/Quantum/DATA/Workbook7F.xlsx']);
T8 = readtable(['~/Dropbox/Directory/Quantum/DATA/Workbook8F.xlsx']);
T9 = readtable(['~/Dropbox/Directory/Quantum/DATA/Workbook9F.xlsx']);
T11 = readtable(['~/Dropbox/Directory/Quantum/DATA/Workbook11F.xlsx']);
T12 = readtable(['~/Dropbox/Directory/Quantum/DATA/Workbook12F.xlsx']);

% Load eigenstates
tptp = readtable('tptp.csv');
tptr1 = readtable('tptr1.csv');
tptr2 = readtable('tptr2.csv');
tr1tp = readtable('tr1tp.csv');
tr1tr1 = readtable('tr1tr1.csv');
tr1tr2 = readtable('tr1tr2.csv');
tr2tp = readtable('tr2tp.csv');
tr2tr1 = readtable('tr2tr1.csv');
tr2tr2 = readtable('tr2tr2.csv');

T2 = T2(1:5,2:6);
T3 = T3(1:5,2:6);
T4 = T4(1:5,2:6);
T5 = T5(1:5,2:6);
T7 = T7(1:5,2:6);
T8 = T8(1:5,2:6);
T9 = T9(1:5,2:6);
T11 = T11(1:5,2:6);
T12 = T12(1:5,2:6);

% Convert to matrix and normalise
T2 = (T2{:,:}.')./sum(T2{:,:},'all');
T3 = T3{:,:}./sum(T3{:,:},'all');
T4 = T4{:,:}./sum(T4{:,:},'all');
T5 = T5{:,:}./sum(T5{:,:},'all');
T7 = T7{:,:}./sum(T7{:,:},'all');
T8 = T8{:,:}./sum(T8{:,:},'all');
T9 = T9{:,:}./sum(T9{:,:},'all');
T11 = T11{:,:}./sum(T11{:,:},'all');
T12 = T12{:,:}./sum(T12{:,:},'all');

T_00 = T2%(T2+T3)./2
T_02 = (T4+T5)./2
T_04 = (T7+T8+T9)./3
T_06 = (T11+T12)./2

I = zeros(3,3);

figure('position',[0,0,1200,500])

range = 2:4
plt = 0
AX = ones(length(range))

for iter = range
    
    eval(['T = T_0' num2str((iter-1)*2)])
        
    % Evaluate overlap integrals
    I(1,1) = sum(tptp{:,:}.*T,'all');
    I(1,2) = sum(tptr1{:,:}.*T,'all');
    I(1,3) = sum(tptr2{:,:}.*T,'all');
    I(2,1) = sum(tr1tp{:,:}.*T,'all');
    I(2,2) = sum(tr1tr1{:,:}.*T,'all');
    I(2,3) = sum(tr1tr2{:,:}.*T,'all');
    I(3,1) = sum(tr2tp{:,:}.*T,'all');
    I(3,2) = sum(tr2tr1{:,:}.*T,'all');
    I(3,3) = sum(tr2tr2{:,:}.*T,'all');


    I

    plt = plt+1
    AX(plt) = subplot(1,length(range),plt)
    b = bar3(I./sum(I, 'all'));
    hold on

    for k = 1:length(b)
        b(k).CData = b(k).ZData;
        b(k).XData = b(k).XData + (min(yvals)-1)*ones(size(b(k).XData));
        b(k).Parent.YDir = 'normal';
    end


    colormap jet; shading interp; axis tight; zlim([0 0.24])
    view([45, 30]);
    xticks([1 2 3]);
    xticklabels({'Tp','Tr1','Tr2'});
    yticks([1 2 3]);
    yticklabels({'Tp','Tr1','Tr2'});
    title(['Disorder = 0.', num2str((iter-1)*2)])
    hold on
    
end