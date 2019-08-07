% loading processed Right Hemisphere human connectome data 
% Roberts method (0.15) on iFOD
load('RH.mat')

%Connectome visualisation
Visualisebrain(RH);


%Degree vs strength plot
binRH = RH>0;
degree = sum(binRH,1)'+sum(binRH,2)-diag(binRH);
degree = degree/2;
strength = sum(RH,1)'+sum(RH,2)-diag(RH);
strengh=strength/2;
subplot(2,2,3);
scatter(degree,strength)
title('B','FontSize', 15);
xlabel('Degree')
ylabel('Strength')

%Degree Histogram
subplot(2,2,4);
hist(degree, 20)
title('C','FontSize', 15);
xlabel('Degree')
ylabel('Frequency')

