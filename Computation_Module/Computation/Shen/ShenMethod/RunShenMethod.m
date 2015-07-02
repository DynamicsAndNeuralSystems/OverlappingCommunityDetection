% ------------------------------------------------------------------------------
% Code for running the Shen et al. method for overlapping community detection
% Relies on maximal clique code, and modularity optimization code.
% Ben Fulcher, 2014-03-20
% Ben Fulcher, 2014-04-17 Implemented sparse calculation, k range, and new
%                           plotting functionality.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% First clear all variables:
% clear all;

% ------------------------------------------------------------------------------
% Set Parameters
% ------------------------------------------------------------------------------

% Set k range to compare: k_r
% k_r = 5; %(6:11);

% Set input matrix
% WhatInput = 'trivial3';
% WhatInput = 'Fig1Shen';
% WhatInput = 'Fig2Shen';
% WhatInput = 'Fig2Shen_weighted';
% WhatInput = 'Brain_1005032';
% WhatInput = 'Brain_1005033';
% WhatInput = 'Brain_1005034';
% WhatInput = 'Alex_F';
WhatInput = 'BenchMark';

% Threshold for counting a link of a weighted network in a resulting binary network
Threshold = 0;

% Whether to save output to a .mat file
DoSaveMat = 1;

% Apply a weight threshold to count maximal cliques for weighted networks
% CliqueWeightThreshold = 1.1; % mean link weight within a clique for it to be counted

% ------------------------------------------------------------------------------
% Define an input matrix:
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
    case 'Fig2Shen_weighted'
        % Weighted version of an example shown in Fig. 2 of the Shen paper:
        A = [0,0.5,1,0.05,0,0,0,0,0,0,0; ...
            0,0,0,0.4,0,0,0,0,0,0,0; ...
            0,0,0,1,0,0.2,1,0,0,0,0; ...
            0,0,0,0,0.6,0,0,0,0,0,0; ...
            0,0,0,0,0,0,0,0,0,0,0; ...
            0,0,0,0,0,0,0,1,0,0,0; ...
            0,0,0,0,0,0,0,0.1,0.3,1,0; ...
            0,0,0,0,0,0,0,0,1,0.9,0.8; ...
            0,0,0,0,0,0,0,0,0,0.2,0; ...
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
        
    case 'Brain_1005032'
        load('1005032_adj_weighted.mat','adj_final')
        A = adj_final;
    case 'Brain_1005033'
        load('1005033_adj_weighted.mat','adj_final')
        A = adj_final;
    case 'Brain_1005034'
        load('1005034_adj_weighted.mat','adj_final')
        A = adj_final;
    case 'Alex_F'
        % Data provided by Alex F, to generate for an online visualization platform
        A = dlmread('Alex_F_Adj.txt');
    case 'BenchMark'
        Ntwk2AdjMat
        A = Adj_Mat;
end

% ------------------------------------------------------------------------------
% Threshold weighted networks
if any(A(A~=0)~=1); % some non-1, non-0 entries (i.e., a weighted network)
    IsWeighted = 1; % a weighted network
    fprintf(1,'ZOMG it''s a weighted network...!\n');
    Aw = A; % keep the weighted version as Aw
    A = (A > Threshold);
    fprintf(1,'After thresholding at %f, we have a link density ~ %f\n',Threshold, ...
        sum(A(:)/length(A)^2));
else
    IsWeighted = 0; % an unweighted network
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
AllMaxClique = maximalCliques(A);
fprintf(1,' Done in %s.\n',BF_thetime(toc(maxcliquetimer)));
clear maxcliquetimer

% AllMaxClique has NumNodes rows, and NumCliques columns
% AllMaxClique contains all maximum cliques, before thresholding with k

% Give some information
fprintf(1,'We found %u maximal cliques in the network\n',size(AllMaxClique,2));


% ------------------------------------------------------------------------------
% Only include maximal cliques with a minimum weight
if IsWeighted && (CliqueWeightThreshold > 0)
    fprintf(1,'Experimenting with applying a clique weight threshold of %f\n', ...
        CliqueWeightThreshold);
    NumCliques = size(AllMaxClique,2);
    MeanCliqueWeight = zeros(NumCliques,1);
    for i = 1:NumCliques
        TheNodes = (AllMaxClique(:,i)>0); % Which nodes are in this clique
        TheLinks = Aw(TheNodes,TheNodes); % Subsection of adjacency matrix
        JustLinks = TheLinks(triu(logical(ones(sum(TheNodes))),+1)); % List of link weights in this clique
        MeanCliqueWeight(i) = mean(JustLinks); % (arithmetic) mean link weight within the clique
    end
    KeepCliques = (MeanCliqueWeight > CliqueWeightThreshold);
    fprintf(1,'Keeping %u / %u maximal cliques with mean link weight exceeding %f!\n', ...
        sum(KeepCliques),NumCliques,CliqueWeightThreshold);
    AllMaxClique = AllMaxClique(:,KeepCliques);
end


% Compute the size of each computed maximal clique
MaxCliqueSize = sum(AllMaxClique);

% ------------------------------------------------------------------------------
% Now we start to use k to threshold the set of all maximal cliques
% Iterate over the full k range specified by the user
% ------------------------------------------------------------------------------

NodeLabels = cell(NumNodes,length(k_r));
QcValues = zeros(length(k_r),1);

for ki = 1:length(k_r)
    
    % Set the k for this run:
    k = k_r(ki);
    
    fprintf(1,'\n------------k = %u-------------\n\n',k);
    
    % Find subordinate maximal cliques
    IsSubordinate = (MaxCliqueSize < k);
    
    fprintf(1,'%u (/%u) maximal clique(s) were large enough to be included as maximal cliques\n', ...
        sum(~IsSubordinate),length(IsSubordinate));
    fprintf(1,['%u (/%u) maximal clique(s) were too small to be included as maximal cliques at k = %u\n' ...
        '----These will be decomposed into subordinate vertices\n'], ...
        sum(IsSubordinate),length(IsSubordinate),k);
    
    % Find nodes that ONLY belong to subordinate maximal cliques
    % (the number of subordinate maximal cliques a vertex belongs to is the same as the number of
    % maximal cliques it belongs to:)
    SubordinateNode = (AllMaxClique*IsSubordinate' == sum(AllMaxClique,2));
    
    fprintf(1,'%u subordinate vertices...\n',sum(SubordinateNode));
    
    % So now we want a network consisting of maximal cliques and subordinate nodes,
    % a cover of the original network
    
    % First make Gprime have only maximal cliques:
    JustMaximalCliques = AllMaxClique(:,~IsSubordinate);
    Gprime = JustMaximalCliques;
    
    % Then append the subordinate vertices:
    if any(SubordinateNode)
        fSubordinateNode = find(SubordinateNode);
        NumSubNode = sum(SubordinateNode);
        % Zeros to start
        AllSubNodes = zeros(NumNodes,NumSubNode);
        % Then add a one in each column for that node:
        for i = 1:NumSubNode
            AllSubNodes(fSubordinateNode(i),i) = 1;
        end
        
        % Add to Gprime
        Gprime = [Gprime,AllSubNodes];
    end
    
    Gprime = logical(Gprime);
    
    % ------------------------------------------------------------------------------
    %% Now we have our Gprime!
    % ------------------------------------------------------------------------------
    CoverSize = size(Gprime,2);
    fprintf(1,'We have %u clusters/subordinate vertices in our cover, Gprime\n',CoverSize);
    
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
    % B = zeros(CoverSize);
    % for i = 1:CoverSize
    %     for j = 1:CoverSize
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
    B = zeros(CoverSize);
    for i = 1:CoverSize
        for j = i:CoverSize
            bigsum = @(x)sum(x(:));
            B(i,j) = bigsum(alpha_s(:,i)*alpha_s(:,j)'.*A_s);
        end
        if CoverSize > 50 && (i==1 || mod(i,floor(CoverSize/10))==0)
            fprintf(1,'Less than %s remaining! We''re at %u / %u\n', ...
                BF_thetime(toc(adjtimer)/i*(CoverSize-i)),i,CoverSize)
        end
    end
    B = B + B'; % Symmetrize
    B(logical(diag(ones(length(B),1)))) = B(logical(diag(ones(length(B),1))))/2; % half diagonal entries
    fprintf(1,' Done in %s.\n',BF_thetime(toc(adjtimer)));
    
    
    % ------------------------------------------------------------------------------
    % Now I just have to run modularity optimization on B
    % ------------------------------------------------------------------------------
    
    modtimer = tic;
    fprintf(1,'Using fast unfolding algorithm for community detection on B...');
    CommunityStruct = cluster_jl(B,1,1,0);
    fprintf(1,' Done in %s :)\n',BF_thetime(toc(modtimer)));
    clear modtimer
    
    
    % ------------------------------------------------------------------------------
    % Extract the communities:
    % ------------------------------------------------------------------------------
    for i = 1:length(CommunityStruct.MOD)
        fprintf(1,'%u: %u communities, modularity = %.3g\n',i, ...
            length(CommunityStruct.SIZE{i}),CommunityStruct.MOD(i));
    end
    
    NumComms = length(CommunityStruct.SIZE{end});
    fprintf(1,'We found %u communities in the cover network after %u iterations\n', ...
        NumComms,length(CommunityStruct.Niter));
    Modularity = CommunityStruct.MOD(end);
    fprintf(1,'With a modularity of %f (min = %f, max = %f)\n', ...
        Modularity, min(CommunityStruct.MOD), max(CommunityStruct.MOD));
    CommLabels = CommunityStruct.COM{end};
    
    % ------------------------------------------------------------------------------
    % Assign community labels to each node in the original graph:
    % ------------------------------------------------------------------------------
    
    % NodeLabels gives each node in the original network a label (or set of labels)
    % indicating its community assignment
    for i = 1:NumNodes
        NodeLabels{i,ki} = unique(CommLabels(Gprime(i,:)));
    end
    
    QcValues(ki) = ComputeQc(ConvertNodeLabelsToCover(NodeLabels(:,ki)),A);
end

% ------------------------------------------------------------------------------
% Save results
% ------------------------------------------------------------------------------
if DoSaveMat
    FileNameSave = sprintf('Shen_%s_k_%u.mat',WhatInput,k);
    save(FileNameSave)
    
    fprintf(1,'Results saved as %s!!\n',FileNameSave);
end

% ------------------------------------------------------------------------------
% Visualize the results using BrainNetViewer
% ------------------------------------------------------------------------------
% PlotCommunities

% ------------------------------------------------------------------------------
% Visualize communities in a simple figure
% ------------------------------------------------------------------------------
% figure('color','w'); box('on');
% for i = 1:length(k_r)
%     subplot(1,length(k_r),i); hold on; box('on')
%     VisualizeCommsBar(NodeLabels(:,i),NumComms,0)
%     % Add a title:
%     NumOverlapping = sum(cellfun(@(x)length(x)>1,NodeLabels(:,i)));
%     title(sprintf(['[k = %u] ' ...
%         '%u comms, %u/%u overlapping, Qc = %.3f.'], ...
%         k_r(i),NumComms,NumOverlapping,NumNodes,QcValues(i)))
%     xlabel('Nodes')
%     ylabel('Community assignments')
% end