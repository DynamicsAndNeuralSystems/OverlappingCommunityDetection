function [NodeMembership] = run_Jerry(AdjMat, NumIter, Thresholds, IsBlind, ListeningRule, NumNodes)

%% Check inputs and set defaults
if nargin < 2 || isempty(NumIter)
    NumIter = 1000;
end
if nargin < 3 || isempty(Thresholds)
    Thresholds = 0.1;
end

%% Code
% Function that runs the Jerry algorithm on a network
% Adapted from code by Ben Fulcher

A = AdjMat;

%% 
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