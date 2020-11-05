close all;
clear;
clc;
%Output=diag([0,1,0,1,0]);
T = readtable('~/Dropbox/Directory/Quantum/DATA/Workbook2F.xlsx');
T = T(1:5,2:6);
T = T{:,:};
load('short-short.mat')
Output = T;
N=length(Output(:,1));
DIM=[N,N];
[V,S,W]=schmidt(Output,DIM);
K=1/sum(S.^4);
sum_S=sqrt(sum(S.^2));
Output=Output/sum_S;
[V2,S2,W2]=schmidt(Output,DIM);
K_norm=1/sum(S2.^4)