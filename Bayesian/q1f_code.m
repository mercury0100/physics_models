clear all
close all
format short

[X,Y] = meshgrid(50:0.1:250,0.05:0.01:0.6);
Posterior = pdf_B(X,Y);
surf(X,Y,Posterior)
drawnow
xlabel X
ylabel Sigma
zlabel f(X,Y)
colormap summer
shading interp
hold on

%Markov Chain starting point
x = 50+rand*200;
y = 0.05+rand*0.45;
z = pdf_B(x,y);

%Set iterations and step size
n = 1000;
deltay = 0.005;
deltax = 2.5;
success=0; %initiate success tally
total=0; %initiate total proposal tally

for i=1:n 
    
    %Proposing a new step 
    nx = x + (randn)*deltax;
    ny = y + (randn)*deltay;
    nz = pdf_B(nx,ny);
    
    %Create the acceptance probability of the proposed step. Since the proposal distribution is normal
    %and therefore symmetric the hastings factor can be reduced to 1.
    ratio = nz/z; %this is the probability ratio between the new proposed/current position
    roll = rand+0.2*ratio > 1; %This generates a 1 with probability min(1,0.2*ratio), otherwise 0
    
    %Saving samples into an array
    samples(i,:) = [nx,ny,nz,ratio];    

    if ratio >= 1
        x = nx;
        y = ny;
        z = nz;
        
        mc(i,:)=[x,y,z];
        success=success+1;
        
        figure(1)
        plot3(nx,ny,nz,'r.')
        drawnow 
        hold on
        
    else
        if roll == 1
            x = nx;
            y = ny;
            z = nz;
            
            mc(i,:)=[x,y,z];
            success=success+1;
            
            figure(1)
            plot3(nx,ny,nz,'r.')
            drawnow 
            hold on
        else
            mc(i,:)=[x,y,z];
        end
    end
   
    total=total+1;
        
end

Acceptance_Rate=success/total

 figure(2)
 h=axes;
 plot(mc(:,3))
 set(h, 'Ydir', 'reverse')
 set(h, 'YAxisLocation', 'Right')
 set(gca,'Ydir','reverse')
 
Mean_X=mean(mc(200:n,1));
StDev_X=std(mc(200:n,1));
Mean_Sigma=mean(mc(200:n,2));
StDev_Sigma=std(mc(200:n,2));
T=table(Acceptance_Rate, Mean_X, StDev_X, Mean_Sigma, StDev_Sigma)

