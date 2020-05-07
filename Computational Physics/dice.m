% dice.m
% Function returning a specified number of random integers, each 
% from one to six.
%
% Inputs:
% ndice - number of random integers required
% Outputs:
% ran6 - row vector with the integers

function ran6 = dice(ndice,nthrows)

ran6=floor(6*rand(nthrows,ndice))+1;
