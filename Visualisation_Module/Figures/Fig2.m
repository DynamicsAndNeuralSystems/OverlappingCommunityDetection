% loading processed Right Hemisphere human connectome data
% Roberts method (0.15) on iFOD
load('RH.mat','RH')

f = figure('color','w');

% Connectome visualisation
subplot(2,3,[1,2,4,5]); axis('square')
Visualisebrain(RH);


%Degree vs strength plot
subplot(2,3,3);
binRH = (RH > 0);
degree = sum(binRH,1)'+sum(binRH,2)-diag(binRH);
degree = degree/2;
strength = sum(RH,1)'+sum(RH,2)-diag(RH);
strength = strength/2;
plot(degree,strength,'ok','MarkerFaceColor','k','MarkerSize',3)
xlim([0,max(degree)])
title('B','FontSize', 15);
xlabel('Degree, k')
ylabel('Strength, s')
% axis('square')

%Degree Histogram
subplot(2,3,6);
axis('square')
histogram(degree,20,'FaceColor',0.5*ones(1,3),'EdgeColor','k')
xlim([0,max(degree)])
title('C','FontSize', 15);
xlabel('Degree, k')
ylabel('Frequency')
% axis('square')

f.Position = [1000        1002         648         336];
