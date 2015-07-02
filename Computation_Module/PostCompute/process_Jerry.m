function [Output] = process_Jerry(NodeMembership, numnodes)
% Function that processes the output of the Jerry method into the matrix

% Finds the number of communities, and makes a matrix of zeros.
Output = zeros(numnodes, max([NodeMembership{:}]));

for i = 1:numnodes
    % equally separates the nodes into communities
    Output(i, NodeMembership{i}) = 1/length(NodeMembership{i});
end