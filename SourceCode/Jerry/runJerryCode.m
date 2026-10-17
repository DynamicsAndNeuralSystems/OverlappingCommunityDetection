function nodeMembership = runJerryCode(adjMatrix,numIter,Thresholds,isBlind,listeningRule,numNodes)
% Overlapping Community Detection (SLPA) algorithm
% isBlind: if true, listeners only sample neighbors' memory from previous iterations;
% if false, they may also use labels already learned earlier in the current sweep.
% Labels below the threshold are dropped, but each node keeps its most frequent label.

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
    nodePosition = zeros(1,numNodes);
    nodePosition(ix) = 1:numNodes; % position of each node in this sweep
    for j = 1:numNodes
        % We focus on this node:
        theNode = ix(j);
        
        % Find neighbors of this node:
        Neighbors = find(adjMatrix(theNode,:) > 0); % the neighbors
        neighborLinkWeights = adjMatrix(theNode,Neighbors); % the link weights to each neighbor
        numNeighbors = length(Neighbors); % the number of neighbors
        
        
        % Isolated nodes have nobody to listen to: keep the previous label
        if numNeighbors==0
            nodeMemory(theNode,i) = nodeMemory(theNode,i-1);
            continue
        end

        % ---LISTENING RULE---
        % Column i of nodeMemory is only filled in for nodes already updated in this
        % sweep (it is still 0 for the rest), so it may only be sampled for those nodes.
        % Blind: only ever sample from previous iterations (columns 1..i-1).
        % Not blind: for neighbors already updated this sweep also include column i.

        switch listeningRule
            case 'majority'
                % Pick the most popular label (weighted)

                % 1) Get a random draw from each neighbor's memory:
                neighborMemorySample = zeros(numNeighbors,1);
                for k = 1:numNeighbors
                    if ~isBlind && nodePosition(Neighbors(k)) < j
                        % Already updated this sweep: its new label (column i) is available
                        neighborMemorySample(k) = nodeMemory(Neighbors(k),randi(i));
                    else
                        % Blind, or not yet updated this sweep: previous iterations only
                        neighborMemorySample(k) = nodeMemory(Neighbors(k),randi(i-1));
                    end
                end

                % 2) Determine the most popular sample from the neighbors, weighting
                % each by the link weights:
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
                if isempty(theNeighbor) % guard against floating-point round-off
                    theNeighbor = Neighbors(end);
                end

                % Pick randomly from that node's memory
                if ~isBlind && nodePosition(theNeighbor) < j
                    % Already updated this sweep: can also pick from the current iteration
                    nodeMemory(theNode,i) = nodeMemory(theNeighbor,randi(i));
                else % Pick only from previous iterations
                    nodeMemory(theNode,i) = nodeMemory(theNeighbor,randi(i-1));
                end
        end

        % A label of 0 (an unfilled memory slot) must never be spoken:
        assert(nodeMemory(theNode,i) > 0,'Invalid (zero) label entered memory')
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
        % As in the original SLPA (Xie, Szymanski & Liu, 2011), labels with frequency
        % below the threshold are removed, but each node keeps at least its most
        % frequent label:
        labelFreq = nodeHistograms{j}(2,:)/max(numIter);
        aboveThreshold = (labelFreq >= Thresholds(i));
        if ~any(aboveThreshold)
            [~,mostFrequent] = max(labelFreq);
            aboveThreshold(mostFrequent) = true;
        end
        nodeMembership{j,i} = nodeHistograms{j}(1,aboveThreshold);
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
