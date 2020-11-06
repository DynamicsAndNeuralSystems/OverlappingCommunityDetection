function alpha = ComputeBelongingCoeff(Cover,adjMat);
% Cover is a matrix with a row for each node, and a column for each community
% A 1 at (i,j) indicates that node i is a member of community j.
%-------------------------------------------------------------------------------
% Ben Fulcher, 2014-03-31
%-------------------------------------------------------------------------------

numNodes = size(Cover,1);
numComms = size(Cover,2);

fprintf(1,'Calculating belonging coefficients for %u nodes in %u communities...', ...
                    numNodes,numComms);

alpha = zeros(numNodes,numComms);

BelongingTimer = tic;

for i = 1:numNodes
    
    if sum(Cover(i,:)==1)
        % Node i is only a member of its own solitary community
        alpha(i,:) = Cover(i,:);
    else
        % Multiple community memberships
        % Need to compute definitions for belonging coefficients...
        
        % Loop over each community to compute belonging coefficients:
        for j = 1:numComms
            % ---Belonging coefficients:
            % How much does node i belong in community j, Eq. (7) in Shen et al., 2009
            % Loop over other nodes in community j
            % if sum(Cover(:,j))==1 % this 'community' in Cover is a subordinate vertex
            %     alpha(i,j) = 0;
            % else
            if (Cover(i,j)==1) % this node is a member of this community
                % Check whether each link that occurs in the full graph, A_ij, occurs in this clique:
                % ij_links_c is a logical expressing this for each link the node has in the full graph,
                % i.e., A(i,:)
                ij_links_c = Cover(adjMat(i,:),j);
                % Check how often the link, A_ij, occurs across all maximal cliques
                % (use all of Cover for now -- the same result, just suboptimal, computationally)
                if sum(ij_links_c) > 0
                    % Looks at columns, corresponding to communities, c, in which node i
                    % is a member (using Cover(~,Cover(i,:)))
                    % Then looks at which of the rows have links to this node in the communities
                    % using Cover(A(i,:)~)
                    % Then takes the sum to count how many communities each link is a member of
                    % [don't really need to look at other communities when ij_links_c is already zero, becase resulting sum will be zero, but hey it's suboptimal...]
                    ij_links_all = sum(Cover(A(i,:),Cover(i,:)),2);
                    % To stop NaNs, we only calculate entries where ij_links_c are nonzero
                    alpha(i,j) = sum(ij_links_c(ij_links_c>0)./ij_links_all(ij_links_c>0));
                % else, alpha(i,j) = 0; % already zero
                end
            % else, alpha(i,j) = 0; % already zero
            end
            % if isnan(alpha(i,j)), keyboard; end
        end
    end
    
    % Normalize entries across this row to sum to 1: [Equivalent to Eq. (8) in Shen et al., 2009]
    % if sum(alpha(i,:)) > 0
    alpha(i,:) = alpha(i,:)/sum(alpha(i,:));

end
fprintf(1,' Done in %s.\n',BF_TheTime(toc(BelongingTimer)));

end
