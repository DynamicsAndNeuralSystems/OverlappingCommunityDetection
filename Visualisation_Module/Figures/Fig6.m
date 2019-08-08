% load RH connectome
load('RH.mat');

%load Louvain community labels for RH connectome
load('louvaincomm.mat');

%Calculate Participation Coefficient
P = participation_coef(RH,louvaincomm);

% Calculate within-module degree z-scored
z = module_degree_zscore(RH, louvaincomm);

scatter(P,z);
hold on

%load OSLOM community affiliation vectors on RH connectome 
load('oslomcomm.mat');
overlapping = find(sum(oslomcomm>0, 2)==2); % list of overlapping nodes obtained from OSLOM

scatter(P(overlapping), z(overlapping), 'r');
xlabel('Participation Coefficient, P');
ylabel('Within-module degree, z');


