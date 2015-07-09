% Ben Fulcher, 2014-04-10
% Unweighted networks

% ------------------------------------------------------------------------------
% First clear all variables:
clear all;

% ------------------------------------------------------------------------------
% Set Parameters
% ------------------------------------------------------------------------------

% ---Set input matrix---

% WhatInput = 'trivial3';
% WhatInput = 'Fig1Shen';
% WhatInput = 'Fig2Shen';
% WhatInput = 'FC';
% WhatInput = 'Brain_1005032';
% WhatInput = 'Brain_1005033';
% WhatInput = 'Brain_1005034';
WhatInput = 'BenchMark';
% WhatInput = 'BenchMarkNorm';

% % ---Number of iterations---
NumIter = 120;

% ---Thresholds---
Thresholds = [0.09];

% ---Algorithmic Parameters---

IsBlind = 1; % Whether node memories update continuously (0), or just after each iterations (1).

% ListeningRule = 'majority'; % Takes a (weighted) majority of labels from neighbors
ListeningRule = 'probabilistic'; % Selects the label from neighbors probabilistically


% ------------------------------------------------------------------------------
% Define an input matrix:
% ------------------------------------------------------------------------------
switch WhatInput
    case 'trivial3'
        A = [0,1,0; 1,0,0; 0,0,0];
        
    case 'Fig2Shen'
        % EXAMPLE SHOWN IN Fig. 2 of Shen paper:
        A = [0,1,1,1,0,0,0,0,0,0,0; ...
            0,0,0,1,0,0,0,0,0,0,0; ...
            0,0,0,1,0,1,1,0,0,0,0; ...
            0,0,0,0,1,0,0,0,0,0,0; ...
            0,0,0,0,0,0,0,0,0,0,0; ...
            0,0,0,0,0,0,0,1,0,0,0; ...
            0,0,0,0,0,0,0,1,1,1,0; ...
            0,0,0,0,0,0,0,0,1,1,1; ...
            0,0,0,0,0,0,0,0,0,1,0; ...
            0,0,0,0,0,0,0,0,0,0,0; ...
            0,0,0,0,0,0,0,0,0,0,0;];
        A = A + A';
    case 'Fig1Shen'
        % EXAMPLE SHOWN IN FIG. 1 OF SHEN PAPER:
        links = [1,2; 1,3; 1,4; 1,5; 1,6; ...
            2,3; 2,4; 2,5; 2,6; ...
            3,4; 3,5; 3,6; 3,7; 3,8; 3,13; ...
            4,5; 4,6; 4,23; ...
            5,6; 5,22; ...
            7,8; 7,9; 7,10; 7,12; 7,13; ...
            8,9; 8,10; 8,11; 8,12; 8,13; ...
            9,10; ...
            10,13; 10,14; 10,15; 10,16; 10,17; ...
            11,12; 11,13; 11,16; 11,17; ...
            12,13; 12,14; 12,16; 12,17; ...
            14,15; 14,16; 14,17; ...
            15,16; 15,17; ...
            16,17; ...
            17,18; 17,19; ...
            18,19; 18,20; 18,21; 18,22; 18,23; ...
            19,22; 19,23; 19,24; ...
            20,21; 20,22; 20,23; ...
            21,23; ...
            22,23; ...
            23,24];
        A = sparse(links(:,1),links(:,2),ones(length(links),1),24,24);
        A = full(A);
        A = A + A';
        
    case 'FC'
        FC = dlmread('weightedFC.txt');
        NumNodes = max(max(FC(:,1)),max(FC(:,2)));
        A = zeros(NumNodes,NumNodes);
        for i = 1:length(FC)
            A(FC(i,1),FC(i,2)) = FC(i,3);
        end
        A = A + A';
        
        
    case 'Brain_1005032'
        load('1005032_adj_weighted.mat','adj_final')
        A = (adj_final>Threshold);
        fprintf(1,'After thresholding at %f, we have a density ~ %f\n',Threshold, ...
            sum(A(:)/length(A)^2));
    case 'Brain_1005033'
        load('1005033_adj_weighted.mat','adj_final')
        A = (adj_final>Threshold);
        fprintf(1,'After thresholding at %f, we have a density ~ %f\n',Threshold, ...
            sum(A(:)/length(A)^2));
    case 'Brain_1005034'
        load('1005034_adj_weighted.mat','adj_final')
        A = (adj_final>Threshold);
        fprintf(1,'After thresholding at %f, we have a density ~ %f\n',Threshold, ...
            sum(A(:)/length(A)^2));
    case 'BenchMark'
        Ntwk2AdjMat;
        A = Adj_Mat;
    case 'BenchMarkNorm'
        Ntwk2AdjMat;
        UnityRescale = @(x) (x-min(x(~isnan(x))))/(max(x(~isnan(x)))-min(x(~isnan(x))));
        % Outlier-adjusted sigmoid:
        Adj_norm = UnityRescale(1./(1 + exp(-(Adj_Mat-median(Adj_Mat(~isnan(Adj_Mat) & Adj_Mat>0))) ...
            /(iqr(Adj_Mat(~isnan(Adj_Mat) & Adj_Mat>0))/1.35))));
        A = Adj_norm;
end


% ------------------------------------------------------------------------------
% A = logical(A);
NumNodes = length(A); % Number of nodes
% ------------------------------------------------------------------------------

fprintf(1,'We have %u nodes, running %u iterations\n', NumNodes,max(NumIter));
if IsBlind
    fprintf(1,'Our listening rule does not incorporate new labels learned during an iteration loop\n');
else
    fprintf(1,'Our listening rule incorporates new labels learned during an iteration loop\n');
end
switch ListeningRule
    case 'majority'
        fprintf(1,'Using a (weighted) majority-of-neighbors listening rule\n');
    case 'probabilistic'
        fprintf(1,'Using a (weighted) probabilistic listening rule\n');
end

AlgorithmTimer = tic;
fprintf(1,'Evaluating speaking/listening rules for %u iterations......\n',NumIter);

% ------------------------------------------------------------------------------
% Run the algorithm
% ------------------------------------------------------------------------------

% Node membership progresses along a row; each row represents that node's memory
NodeMemory = zeros(NumNodes,max(NumIter));

% 1. Assign initial node memberships
NodeMemory(:,1) = (1:NumNodes)';

for i = 2:max(NumIter)
    % Determine the permutation to shuffle through nodes:
    ix = randperm(NumNodes);
    for j = 1:NumNodes
        % We focus on this node:
        TheNode = ix(j);
        
        % Find neighbors of this node:
        Neighbors = find(A(TheNode,:)>0); % the neighbors
        NeighborLinkWeights = A(TheNode,Neighbors); % the link weights to each neighbor
        NumNeighbors = length(Neighbors); % the number of neighbors
        
        
        % ---LISTENING RULE---
        
        switch ListeningRule
            case 'majority'
                % Pick the most popular label (weighted)
                
                % 1) Get a random draw from each neighbor's memory:
                NeighborMemorySample = zeros(NumNeighbors,1);
                if IsBlind
                    % Compute samples of neighbors based on previous iterations
                    % (even if a node has already received a new label during this iteration)
                    r = randi(i-1,[NumNeighbors,1]); % Pick a random element from each neighbor
                    for k = 1:NumNeighbors
                        NeighborMemorySample(k) = NodeMemory(Neighbors(k),randi(i-1)); % NodeMemory(k,r(k));
                    end
                else
                    for k = 1:NumNeighbors
                        if find(ix==Neighbors(k),1) < j
                            NeighborMemorySample(k) = NodeMemory(Neighbors(k),randi(i)); % NodeMemory(k,r(k));
                        else % already selected this turn - use the new part of memory
                            NeighborMemorySample(k) = NodeMemory(Neighbors(k),randi(i-1)); % NodeMemory(k,r(k));
                        end
                    end
                end
                
                % 2) Pick a number from this sample, by weighting each by the link weights:
                % Determine the most popular sample from the neighbors
                UniqueLabels = unique(NeighborMemorySample);
                LabelWeights = zeros(length(UniqueLabels),1);
                for k = 1:length(UniqueLabels)
                    LabelWeights(k) = sum(NeighborLinkWeights(NeighborMemorySample==UniqueLabels(k)));
                end
                
                % 3) Pick the most popular (weighted) label:
                % If multiple values have maximum weight, select one at random
                MostPopularLabel = UniqueLabels(LabelWeights==max(LabelWeights));
                if length(MostPopularLabel)==1 % Just one popular one: assign it
                    NodeMemory(TheNode,i) = MostPopularLabel;
                else
                    % Select one (from the most popular labels) at random
                    NodeMemory(TheNode,i) = MostPopularLabel(randi(length(MostPopularLabel)));
                end
                
            case 'probabilistic'
                % Get a random draw from the memory of a randomly selected (weighted)
                % neighbor's memory
                
                % First, pick a random neighbor:
                r = rand*sum(NeighborLinkWeights);
                TheNeighbor = Neighbors(find(r <= cumsum(NeighborLinkWeights),1,'first'));
                
                % Pick randomly from that node's memory
                if IsBlind % Potentially pick from current iteration
                    NodeMemory(TheNode,i) = NodeMemory(TheNeighbor,randi(i));
                else % Pick only from previous iterations
                    NodeMemory(TheNode,i) = NodeMemory(TheNeighbor,randi(i-1));
                end
        end
    end
end

fprintf(1,'Speaking-listening rules completed for %u iterations in %s\n', ...
    NumIter,BF_thetime(toc(AlgorithmTimer)));
clear AlgorithmTimer

% ------------------------------------------------------------------------------
% Compute histograms for each node's memory
% ------------------------------------------------------------------------------

NumThresholds = length(Thresholds);

% First generate histograms for each node's memory using the 'unique' function

NodeHistograms = cell(NumNodes,1);
% Each element is a set of labels in the node's memory, and its frequency.
for i = 1:NumNodes
    LabelsInMemory = unique(NodeMemory(i,:));
    NodeHistograms{i} = zeros(2,length(LabelsInMemory));
    NodeHistograms{i}(1,:) = LabelsInMemory;
    NodeHistograms{i}(2,:) = arrayfun(@(x)sum(NodeMemory(i,:)==x),LabelsInMemory);
end

% ------------------------------------------------------------------------------
% Assign community memberships using thresholds
% ------------------------------------------------------------------------------

NodeMembership = cell(NumNodes,NumThresholds);
NumComms = zeros(NumThresholds,1);

for i = 1:NumThresholds
    for j = 1:NumNodes
        AboveThreshold = NodeHistograms{j}(2,:)/NumIter >= Thresholds(i);
        if sum(AboveThreshold) > 0
            NodeMembership{j,i} = NodeHistograms{j}(1,AboveThreshold);
        else
            NodeMembership{j,i} = [];
        end
    end
    AllCommunityLabels = horzcat(NodeMembership{:,i});
    
    % Now label each community according to the biggest (label 1), 2, 3, ... etc.
    % down to the smallest.
    UniqueCommLabs = unique(AllCommunityLabels);
    LabelFrequency = arrayfun(@(x)sum(AllCommunityLabels==x),UniqueCommLabs);
    [~,ix] = sort(LabelFrequency,'descend');
    SortedCommLabels = UniqueCommLabs(ix); % these are the sorted community labels
    
    % Now we want to loop through and rename all nodes with the new community labels
    for j = 1:NumNodes
        NodeMembership{j,i} = arrayfun(@(x)find(SortedCommLabels==x),NodeMembership{j,i});
    end
    
    NumOverlapping = sum(cellfun(@(x)length(x)>1,NodeMembership(:,i)));
    NumUnassigned = sum(cellfun(@(x)length(x)==0,NodeMembership(:,i)));
    
    NumComms(i) = length(UniqueCommLabs);
    
    fprintf(1,['At a threshold %f (%.1f/%u), we have %u communities, ' ...
        '%u overlapping nodes, %u nodes without a community.\n'], ...
        Thresholds(i),Thresholds(i)*NumIter,NumIter,NumComms(i),NumOverlapping,NumUnassigned);
end


% ------------------------------------------------------------------------------
% Visualize
% ------------------------------------------------------------------------------
BarHeight = 1;
BarWidth = 1;

figure('color','w'); box('on');
for i = 1:NumThresholds
    subplot(NumThresholds,1,i); box('on')
    title(sprintf('r = %f---%u communities',Thresholds(i),NumComms(i)))
    
    
    if NumComms(i) <= 9
        Colors = BF_getcmap('set1',NumComms(i),0);
    elseif NumComms(i) <= 21
        Colors = [BF_getcmap('set1',NumComms(i),0);BF_getcmap('set3',NumComms(i),0)];
    elseif NumComms(i) <= 64
        Colors = jet(NumComms(i));
    else
        fprintf('%u is too many communities for me to color!\n',NumComms(i));
        continue
    end
    
    for j = 1:NumComms(i)
        for k = 1:NumNodes
            if any(NodeMembership{k,i}==j);
                rectangle('Position',[k-BarWidth/2,j-BarHeight/2,BarWidth,BarHeight],'FaceColor', ...
                    Colors(j,:),'EdgeColor',Colors(j,:),'LineWidth',0.01)
            end
        end
    end
    xlim([0.5,NumNodes+0.5])
    if NumComms(i)>0
        ylim([0.5,NumComms(i)+0.5])
    end
end