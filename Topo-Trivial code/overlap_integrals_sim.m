
% Import experimental data table
T = readtable('sim_output_0.8.csv');
T = T(1:5,1:5);
% Convert to matrix and normalise
T = T{:,:};
T = T./sum(T,'all');

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

% Evaluate overlap integrals
I = zeros(3,3);

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

b = bar3(I./sum(I,'all'));

%for k = 1:length(b)
%    b(k).CData = b(k).ZData;
%    b(k).XData = b(k).XData + (min(yvals)-1)*ones(size(b(k).XData));
%    b(k).Parent.YDir = 'normal';
%end

colormap jet; shading interp; colorbar; axis tight;
view([45, 30]);
xticks([1 2 3]);
xticklabels({'Tp','Tr1','Tr2'});
yticks([1 2 3]);
yticklabels({'Tp','Tr1','Tr2'});