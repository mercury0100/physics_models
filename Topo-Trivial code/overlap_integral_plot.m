
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

levels = [0, 0.2, 0.4, 0.6, 0.8];
I = zeros(3,3, length(levels));

for i = levels
    % Import experimental data table
    str = ['avg_sim_output_', num2str(i), '.csv']
    T = readtable(str);
    T = T(1:5,1:5);
    % Convert to matrix and normalise
    T = T{:,:}.^2;
    T = T./sum(T,'all');

    % Evaluate overlap integrals
    idx = i*5 + 1

    I(1,1,idx) = sum(tptp{:,:}.*T,'all');
    I(1,2,idx) = sum(tptr1{:,:}.*T,'all');
    I(1,3,idx) = sum(tptr2{:,:}.*T,'all');
    I(2,1,idx) = sum(tr1tp{:,:}.*T,'all');
    I(2,2,idx) = sum(tr1tr1{:,:}.*T,'all');
    I(2,3,idx) = sum(tr1tr2{:,:}.*T,'all');
    I(3,1,idx) = sum(tr2tp{:,:}.*T,'all');
    I(3,2,idx) = sum(tr2tr1{:,:}.*T,'all');
    I(3,3,idx) = sum(tr2tr2{:,:}.*T,'all');
    
    I(:,:,idx) = I(:,:,idx)./sum(I(:,:,idx),'all');

    I;
end

tp_tp(:) = I(1,1,:);
tp_tr1(:) = I(1,2,:);
tp_tr2(:) = I(1,3,:);
tr1_tp(:) = I(2,1,:);
tr1_tr1(:) = I(2,2,:);
tr1_tr2(:) = I(2,3,:);
tr2_tp(:) = I(3,1,:);
tr2_tr1(:) = I(3,2,:);
tr2_tr2(:) = I(3,3,:);

figure(1)
hold on
plot(levels, tp_tp, 'b--o');
plot(levels, tp_tr1);
plot(levels, tp_tr2);
plot(levels, tr1_tp);
plot(levels, tr1_tr1, '--x');
plot(levels, tr1_tr2);
plot(levels, tr2_tp);
plot(levels, tr2_tr1);
plot(levels, tr2_tr2);
legend('tptp', 'tptr1', 'tptr2', 'tr1tp', 'tr1tr1', 'tr1tr2', 'tr2tp', 'tr2tr1', 'tr2tr2');
title('Overlap integrals of Average profile with entangled eigenstates vs Disorder');