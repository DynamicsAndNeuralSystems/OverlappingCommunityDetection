function nodeLabels = runShenCode(AdjMat, cSize, cwThresh, lThresh)
% Function that runs the Shen method on the matrix
% Adapted from code by Ben Fulcher
%-------------------------------------------------------------------------------

A = AdjMat; % Adjacency matrix
cliqueWeightThreshold = cwThresh; 
Threshold = lThresh; % Link threshold
k_r = cSize; % Clique size

% ------------------------------------------------------------------------------
% Threshold weighted networks
if any(A(A~=0)~=1); % some non-1, non-0 entries (i.e., a weighted network)
    isWeighted = true; % a weighted network
    fprintf(1,'ZOMG it''s a weighted network...!\n');
    Aw = A; % keep the weighted version as Aw
    A = (A > Threshold);
    fprintf(1,'After thresholding at %f, we have a link density ~ %f\n',Threshold, ...
        sum(A(:)/length(A)^2));
else
    isWeighted = false; % an unweighted network
end

A = logical(A); % binary, logical matrix
A_s = sparse(A); % make a sparse version, can speed up some calculations later
NumNodes = length(A);

% ------------------------------------------------------------------------------
fprintf(1,'Finding overlapping communities using the method of Shen et al. for k = ');
for i = 1:length(k_r)
    if i < length(k_r)
        fprintf(1,'%u, ',k_r(i));
    else
        fprintf(1,'%u.\n',k_r(i));
    end
end

% ------------------------------------------------------------------------------
% Find the maximal cliques in the input matrix:
maxcliquetimer = tic;
fprintf(1,'Finding maximal cliques in a %ux%u network using maximalCliques code...',length(A),length(A));

allMaxClique = maximalCliques(A);
fprintf(1,' Done in %s.\n',BF_TheTime(toc(maxcliquetimer)));
clear maxcliquetimer

% allMaxClique has NumNodes rows, and numCliques columns
% allMaxClique contains all maximum cliques, before thresholding with k

% Give some information
fprintf(1,'We found %u maximal cliques in the network\n',size(allMaxClique,2));


% ------------------------------------------------------------------------------
% Only include maximal cliques with a minimum weight
if isWeighted && (cliqueWeightThreshold > 0)
    fprintf(1,'Experimenting with applying a clique weight threshold of %f\n', ...
        cliqueWeightThreshold);
    numCliques = size(allMaxClique,2);
    meanCliqueWeight = zeros(numCliques,1);
    for i = 1:numCliques
        theNodes = (allMaxClique(:,i)>0); % Which nodes are in this clique
        theLinks = Aw(theNodes,theNodes); % Subsection of adjacency matrix
        justLinks = theLinks(triu(logical(ones(sum(theNodes))),+1)); % List of link weights in this clique
        meanCliqueWeight(i) = mean(justLinks); % (arithmetic) mean link weight within the clique
    end
    keepCliques = (meanCliqueWeight > cliqueWeightThreshold);
    fprintf(1,'Keeping %u / %u maximal cliques with mean link weight exceeding %f!\n', ...
        sum(keepCliques),numCliques,cliqueWeightThreshold);
    allMaxClique = allMaxClique(:,keepCliques);
end


% Compute the size of each computed maximal clique
maxCliqueSize = sum(allMaxClique);

% ------------------------------------------------------------------------------
% Now we start to use k to threshold the set of all maximal cliques
% Iterate over the full k range specified by the user
% ------------------------------------------------------------------------------

nodeLabels = cell(NumNodes,length(k_r));
QcValues = zeros(length(k_r),1);

for ki = 1:length(k_r)
    
    % Set the k for this run:
    k = k_r(ki);
    
    fprintf(1,'\n------------k = %u-------------\n\n',k);
    
    % Find subordinate maximal cliques
    isSubordinate = (maxCliqueSize < k);
    
    fprintf(1,'%u (/%u) maximal clique(s) were large enough to be included as maximal cliques\n', ...
        sum(~isSubordinate),length(isSubordinate));
    fprintf(1,['%u (/%u) maximal clique(s) were too small to be included as maximal cliques at k = %u\n' ...
        '----These will be decomposed into subordinate vertices\n'], ...
        sum(isSubordinate),length(isSubordinate),k);
    
    % Find nodes that ONLY belong to subordinate maximal cliques
    % (the number of subordinate maximal cliques a vertex belongs to is the same as the number of
    % maximal cliques it belongs to:)
    subordinateNode = (allMaxClique*isSubordinate' == sum(allMaxClique,2));
    
    fprintf(1,'%u subordinate vertices...\n',sum(subordinateNode));
    
    % So now we want a network consisting of maximal cliques and subordinate nodes,
    % a cover of the original network
    
    % First make Gprime have only maximal cliques:
    justMaximalCliques = allMaxClique(:,~isSubordinate);
    Gprime = justMaximalCliques;
    
    % Then append the subordinate vertices:
    if any(subordinateNode)
        fsubordinateNode = find(subordinateNode);
        numSubNode = sum(subordinateNode);
        % Zeros to start
        allSubNodes = zeros(NumNodes,numSubNode);
        % Then add a one in each column for that node:
        for i = 1:numSubNode
            allSubNodes(fsubordinateNode(i),i) = 1;
        end
        
        % Add to Gprime
        Gprime = [Gprime,allSubNodes];
    end
    
    Gprime = logical(Gprime);
    
    % ------------------------------------------------------------------------------
    %% Now we have our Gprime!
    % ------------------------------------------------------------------------------
    coverSize = size(Gprime,2);
    fprintf(1,'We have %u clusters/subordinate vertices in our cover, Gprime\n',coverSize);
    
    % ------------------------------------------------------------------------------
    %% Compute Belonging Coefficients for the cover specified by Gprime
    % ------------------------------------------------------------------------------
    
    alpha = ComputeBelongingCoeff(Gprime,A);
    
    % ------------------------------------------------------------------------------
    % So we have belonging coefficients
    % Now we need to construct our new adjacency matrix for the cover of the network
    % This is Eq. (9) in Shen et al. J. Stat. Mech. (2009)
    % ------------------------------------------------------------------------------
    adjtimer = tic;
    
    % ---Implementation 1---
    % fprintf(1,'Constructing new adjacency matrix from belonging coefficients...');
    % B = zeros(coverSize);
    % for i = 1:coverSize
    %     for j = 1:coverSize
    %         bigsum = @(x)sum(x(:));
    %         B(i,j) = bigsum(alpha(:,i)*alpha(:,j)'.*A);
    %     end
    % end
    % fprintf(1,' Done in %s.\n',toc(adjtimer));
    
    % ---Implementation 2---
    % This version should be a little more efficient for large, sparse networks
    % Just goes through half the links then symmetrizes...
    % Sparse version is about twice as fast
    alpha_s = sparse(alpha);
    fprintf(1,'Constructing a new adjacency matrix from belonging coefficients...');
    B = zeros(coverSize);
    for i = 1:coverSize
        for j = i:coverSize
            bigsum = @(x)sum(x(:));
            B(i,j) = bigsum(alpha_s(:,i)*alpha_s(:,j)'.*A_s);
        end
        if coverSize > 50 && (i==1 || mod(i,floor(coverSize/10))==0)
            fprintf(1,'Less than %s remaining! We''re at %u / %u\n', ...
                BF_TheTime(toc(adjtimer)/i*(coverSize-i)),i,coverSize)
        end
    end
    B = B + B'; % Symmetrize
    B(logical(diag(ones(length(B),1)))) = B(logical(diag(ones(length(B),1))))/2; % half diagonal entries
    fprintf(1,' Done in %s.\n',BF_TheTime(toc(adjtimer)));
    
    
    % ------------------------------------------------------------------------------
    % Now I just have to run modularity optimization on B
    % ------------------------------------------------------------------------------
    
    modtimer = tic;
    fprintf(1,'Using fast unfolding algorithm for community detection on B...');
    communityStruct = cluster_jl(B,1,1,0);
    fprintf(1,' Done in %s :)\n',BF_TheTime(toc(modtimer)));
    clear modtimer
    
    
    % ------------------------------------------------------------------------------
    % Extract the communities:
    % ------------------------------------------------------------------------------
    for i = 1:length(communityStruct.MOD)
        fprintf(1,'%u: %u communities, modularity = %.3g\n',i, ...
            length(communityStruct.SIZE{i}),communityStruct.MOD(i));
    end
    
    NumComms = length(communityStruct.SIZE{end});
    fprintf(1,'We found %u communities in the cover network after %u iterations\n', ...
        NumComms,length(communityStruct.Niter));
    Modularity = communityStruct.MOD(end);
    fprintf(1,'With a modularity of %f (min = %f, max = %f)\n', ...
        Modularity, min(communityStruct.MOD), max(communityStruct.MOD));
    CommLabels = communityStruct.COM{end};
    
    % ------------------------------------------------------------------------------
    % Assign community labels to each node in the original graph:
    % ------------------------------------------------------------------------------
    
    % nodeLabels gives each node in the original network a label (or set of labels)
    % indicating its community assignment
    for i = 1:NumNodes
        nodeLabels{i,ki} = unique(CommLabels(Gprime(i,:)));
    end
    
    QcValues(ki) = ComputeQc(ConvertNodeLabelsToCover(nodeLabels(:,ki)),A);
end
