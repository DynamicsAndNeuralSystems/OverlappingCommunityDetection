% load RH connectome
% load('RH.mat');

%load Louvain community labels for RH connectome
load('louvaincomm.mat');

%Calculate Participation Coefficient
P = participation_coef(RH,louvaincomm);

% Calculate within-module degree z-scored
z = module_degree_zscore(RH, louvaincomm);

% subplot(1,2,1); axis square;

color1 = [0.66,0.66,0.66];
scatter(P,z,80,color1,'filled');
hold on

%load OSLOM community affiliation vectors on RH connectome 
load('oslomcomm.mat');
overlapping = find(sum(oslomcomm>0, 2)==2); % list of overlapping nodes obtained from OSLOM

scatter(P(overlapping), z(overlapping), 80,'filled','r');
xlabel('Participation Coefficient, P');
ylabel('Within-module strength, z');
xlim([0 1]);
legend('Non-Overlapping nodes','Overlapping nodes','Location','southeast');
title('Louvain', 'Fontsize',15);

% For OSLOM, replacing each overlapping node by two non-overlapping nodes
% (since every overlapping node belong to two communities)

RHnew(1:180,1:180)=RH;
RHnew(181:191,1:180)=RH(overlapping,:);
RHnew(1:180,181:191)=RH(:,overlapping);
RHnew(181:191,181:191) = RH(overlapping, overlapping);

%load OSLOM community labels for RHnew connectome
load('oslomRHnew.mat');

%Calculate Participation Coefficient
P = participation_coef(RHnew,oslomRHnew);

% Calculate within-module degree z-scored
z = module_degree_zscore(RHnew, oslomRHnew);

% subplot(1,2,2); axis square;
figure;
color1 = [0.66,0.66,0.66];
scatter(P,z,80,color1,'filled');
hold on

overlapping(12:22) = [181:1:191];

scatter(P(overlapping), z(overlapping), 80,'filled','r');
xlabel('Participation Coefficient, P');
ylabel('Within-module strength, z');
xlim([0 1]);
legend('Non-Overlapping nodes','Overlapping nodes','Location','southeast');
title('OSLOM', 'Fontsize',15);




