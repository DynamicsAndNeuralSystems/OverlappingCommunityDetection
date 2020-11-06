function nodeMembership = runJerryCode(adjMatrix,numIter,Thresholds,isBlind,listeningRule,numNodes)
% Overlapping Community Detection (SLPA) algorithm
% 

%-------------------------------------------------------------------------------
%% Check inputs and set defaults
%-------------------------------------------------------------------------------
if nargin < 2 || isempty(numIter)
    numIter = 1000;
end
if nargin < 3 || isempty(Thresholds)
    Thresholds = 0.1;
end
if nargin < 4
    isBlind = true;
end
if nargin < 5
    listeningRule = 'probabilistic';
end
if nargin < 6
    numNodes = size(adjMatrix,1);
end
%-------------------------------------------------------------------------------

%% 
fprintf(1,'We have %u nodes, running %u iterations\n',numNodes,max(numIter));
if isBlind
    fprintf(1,'Our listening rule does not incorporate new labels learned during an iteration loop\n');
else
    fprintf(1,'Our listening rule incorporates new labels learned during an iteration loop\n');
end
switch listeningRule
    case 'majority'
        fprintf(1,'Using a (weighted) majority-of-neighbors listening rule\n');
    case 'probabilistic'
        fprintf(1,'Using a (weighted) probabilistic listening rule\n');
    otherwise
        error('Unknown listening rule, %s',listeningRule)
end

algorithmTimer = tic;
fprintf(1,'Evaluating speaking/listening rules for %u iterations......\n',numIter);

% ------------------------------------------------------------------------------
% Run the algorithm
% ------------------------------------------------------------------------------

% Node membership progresses along a row; each row represents that node's memory
nodeMemory = zeros(numNodes,max(numIter));

% 1. Assign initial node memberships
nodeMemory(:,1) = (1:numNodes)';

for i = 2:max(numIter)
    % Determine the permutation to shuffle through nodes:
    ix = randperm(numNodes);
    for j = 1:numNodes
        % We focus on this node:
        theNode = ix(j);
        
        % Find neighbors of this node:
        Neighbors = find(adjMatrix(theNode,:) > 0); % the neighbors
        neighborLinkWeights = adjMatrix(theNode,Neighbors); % the link weights to each neighbor
        numNeighbors = length(Neighbors); % the number of neighbors
        
        
        % ---LISTENING RULE---
        
        switch listeningRule
            case 'majority'
                % Pick the most popular label (weighted)
                
                % 1) Get a random draw from each neighbor's memory:
                neighborMemorySample = zeros(numNeighbors,1);
                if isBlind
                    % Compute samples of neighbors based on previous iterations
                    % (even if a node has already received a new label during this iteration)
                    r = randi(i-1,[numNeighbors,1]); % Pick a random element from each neighbor
                    for k = 1:numNeighbors
                        neighborMemorySample(k) = nodeMemory(Neighbors(k),randi(i-1)); % nodeMemory(k,r(k));
                    end
                else
                    for k = 1:numNeighbors
                        if find(ix==Neighbors(k),1) < j
                            neighborMemorySample(k) = nodeMemory(Neighbors(k),randi(i)); % nodeMemory(k,r(k));
                        else % already selected this turn - use the new part of memory
                            neighborMemorySample(k) = nodeMemory(Neighbors(k),randi(i-1)); % nodeMemory(k,r(k));
                        end
                    end
                end
                
                % 2) Pick a number from this sample, by weighting each by the link weights:
                % Determine the most popular sample from the neighbors
                uniqueLabels = unique(neighborMemorySample);
                labelWeights = zeros(length(uniqueLabels),1);
                for k = 1:length(uniqueLabels)
                    labelWeights(k) = sum(neighborLinkWeights(neighborMemorySample==uniqueLabels(k)));
                end
                
                % 3) Pick the most popular (weighted) label:
                % If multiple values have maximum weight, select one at random
                mostPopularLabel = uniqueLabels(labelWeights==max(labelWeights));
                if length(mostPopularLabel)==1 % Just one popular one: assign it
                    nodeMemory(theNode,i) = mostPopularLabel;
                else
                    % Select one (from the most popular labels) at random
                    nodeMemory(theNode,i) = mostPopularLabel(randi(length(mostPopularLabel)));
                end
                
            case 'probabilistic'
                % Get a random draw from the memory of a randomly selected (weighted)
                % neighbor's memory
                
                % First, pick a random neighbor:
                r = rand*sum(neighborLinkWeights);
                theNeighbor = Neighbors(find(r <= cumsum(neighborLinkWeights),1,'first'));
                
                % Pick randomly from that node's memory
                if isBlind % Potentially pick from current iteration
                    nodeMemory(theNode,i) = nodeMemory(theNeighbor,randi(i));
                else % Pick only from previous iterations
                    nodeMemory(theNode,i) = nodeMemory(theNeighbor,randi(i-1));
                end
        end
    end
end

fprintf(1,'Speaking-listening rules completed for %u iterations in %s\n', ...
        numIter,BF_TheTime(toc(algorithmTimer)));
clear('algorithmTimer')

% ------------------------------------------------------------------------------
% Compute histograms for each node's memory
% ------------------------------------------------------------------------------

numThresholds = length(Thresholds);

% First generate histograms for each node's memory using the 'unique' function

nodeHistograms = cell(numNodes,1);
% Each element is a set of labels in the node's memory, and its frequency.
for i = 1:numNodes
    labelsInMemory = unique(nodeMemory(i,:));
    nodeHistograms{i} = zeros(2,length(labelsInMemory));
    nodeHistograms{i}(1,:) = labelsInMemory;
    nodeHistograms{i}(2,:) = arrayfun(@(x)sum(nodeMemory(i,:)==x),labelsInMemory);
end

% ------------------------------------------------------------------------------
% Assign community memberships using thresholds
% ------------------------------------------------------------------------------

nodeMembership = cell(numNodes,numThresholds);
numComms = zeros(numThresholds,1);

for i = 1:numThresholds
    for j = 1:numNodes
        aboveThreshold = nodeHistograms{j}(2,:)/numIter >= Thresholds(i);
        if sum(aboveThreshold) > 0
            nodeMembership{j,i} = nodeHistograms{j}(1,aboveThreshold);
        else
            nodeMembership{j,i} = [];
        end
    end
    allCommunityLabels = horzcat(nodeMembership{:,i});
    
    % Now label each community according to the biggest (label 1), 2, 3, ... etc.
    % down to the smallest.
    uniqueCommLabs = unique(allCommunityLabels);
    labelFrequency = arrayfun(@(x)sum(allCommunityLabels==x),uniqueCommLabs);
    [~,ix] = sort(labelFrequency,'descend');
    lortedCommLabels = uniqueCommLabs(ix); % these are the sorted community labels
    
    % Now we want to loop through and rename all nodes with the new community labels
    for j = 1:numNodes
        nodeMembership{j,i} = arrayfun(@(x)find(lortedCommLabels==x),nodeMembership{j,i});
    end
    
    numOverlapping = sum(cellfun(@(x)length(x)>1,nodeMembership(:,i)));
    numUnassigned = sum(cellfun(@(x)length(x)==0,nodeMembership(:,i)));
    
    numComms(i) = length(uniqueCommLabs);
    
    fprintf(1,['At a threshold %f (%.1f/%u), we have %u communities, ' ...
        '%u overlapping nodes, %u nodes without a community.\n'], ...
        Thresholds(i),Thresholds(i)*numIter,numIter,numComms(i),numOverlapping,numUnassigned);
end
