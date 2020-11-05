%clear all
%load('/Users/ablanco/Desktop/Schmidt/topo40.mat')
% load('/Users/ablanco/Desktop/Schmidt/topo40plotted.mat')
%load('/Users/ablanco/Desktop/Schmidt/topo60.mat')
% load('/Users/ablanco/Desktop/Schmidt/topo60plotted.mat')
% load('/Users/ablanco/Desktop/Schmidt/triv00.mat')
load('/Users/ablanco/Desktop/Schmidt/triv00plotted.mat')
% load('/Users/ablanco/Desktop/Schmidt/triv40.mat')
%load('/Users/ablanco/Desktop/Schmidt/triv40plotted.mat')
x = sqrt(H);
y = svd(x);
z = sum(y.^4);
K = (1/z)*(sum(y.^2))^2
