close all 
clear all
format long

% Defining the Data vectors d and z and the known deviation sigma
Data = [0.19,12;0.28,29;0.37,31;0.46,45;0.55,53;0.64,50;0.73,46;0.82,54;0.91,59;1.00,83;1.09,78;1.18,88;1.27,104;1.36,112;1.45,124;1.54,146;1.63,128;1.72,128;1.81,181;1.90,150;1.99,169;2.08,162;2.17,189;2.26,183;2.34,201;2.43,201;2.52,219;2.61,229;2.70,233;2.79,248;2.88,255;2.97,232;3.06,262;3.15,279;3.24,281;3.33,300;3.42,311;3.51,309;3.60,331;3.69,328;3.78,335;3.87,361;3.96,382;4.05,368;4.14,385;4.23,428;4.32,422;4.41,429;4.50,418];
sigma = 25;     %(km/s)
d = Data(:,1);  %(km/s)
z = Data(:,2);  %(Mpc)

%Initialize active sample of size N, Evidence and X(i)
N = 100;
Evidence = 0;
X(1) = 1;
active = zeros(N,2);

%Generates an active set of N likelihood values with corresponding A values
%from the logarithmic prior
for n=1:N 
    active(n,:) = likelihood_2(sigma,d,z);
end

for i=2:1001
    
    %Select the lowest likelihood point and integrate
    active = sortrows(active,1);
    LogL = active(1,1);
    A = active(1,2);
    L = exp(active(1,1));
    X(i) = exp(-(i-1)/N);
    width = X(i-1) - X(i);
    Delta_E = L*width;
    Evidence = Evidence + Delta_E;
    LogLnew=LogL-1;
    
    %Replace the lowest likelihood point with a random selection from the p
    %rior volume with higher likelihood
    while LogLnew < LogL
       active(1,:) = likelihood_2(sigma,d,z);
       LogLnew = active(1,1);
    end
    
end 

LogEvidence=log(Evidence);
T = table(A, Evidence, LogEvidence)

%Plot the redshift model from optimised A
x = [0:0.1:5];
y = A*x;

%Plot linear fit and data with error bars
figure(1)
err = sigma*ones(length(d),1).';
errorbar(d,z,err)
hold on
plot(x,y)

%Plot residuals with error bars
figure(2)
residuals=z-A*d;
errorbar(residuals,err,'bs')