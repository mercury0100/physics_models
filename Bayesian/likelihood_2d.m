function Output = likelihood_2d(sigma,d,z)

    %Number of points
    DataPoints = length(d);
    Ao = z(:,1)./d(:,1);
    
    %This finds the range of A 
    LogMinA = log(min(Ao));
    LogMaxA = log(max(Ao));
    
    %Select a random A in log space and a random A2 from a uniform prior
    %on [0,10];
    LogA = LogMinA+(LogMaxA-LogMinA)*rand;
    A = exp(LogA);
    A2 = rand*10;
    sum = 0;
    
    %Generate the probability of A given the data and gaussian error with
    %known deviation
    for i=1:DataPoints
        ExponentialTerm = ((A*d(i)+A2*d(i)^2 - z(i))^2)/(2*sigma^2);
        sum = sum + ExponentialTerm;
    end
    FrontFactor = DataPoints*log(sqrt(2*pi)*sigma);
    LogL = -sum-FrontFactor;
    Output = [LogL,A,A2];
end