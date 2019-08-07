% load ENMI for each method for 1000 benchmarks
clear
cd 'ENMI'
load('ENMI_OSLOM.mat');
load('ENMI_NNMF.mat'); 
load('ENMI_Shen.mat');  
load('ENMI_Infomap.mat'); 
load('ENMI_Jerry.mat'); 
cd ../
%-------------------------------------------------------------------------

% Sort performance of methods by mean
methods = [oslom, shen, nnmf, info, jerry(:,5)]; %only using jerry at 0.5 threshold as it is the suggested value and ENMI values across all thresholds are fairly similar
meanarray = mean(methods);
[meanarray, idx] = sort(-meanarray);
methods = methods(:,idx);
[m,n] = size(methods);
methods = mat2cell(methods,m,ones(n,1));

names = {'OSLOM\_50','Clique\_3','Clique\_4','Clique\_5','Clique\_6','Clique\_7','Clique\_9','NNMF\_10','NNMF\_20','NNMF\_30','NNMF\_40','Infomap','SLPA'};
names = names(idx);


% Call function to create violin plots
BF_JitteredParallelScatter(methods,names);