function [Output] = process_Shen(NodeLabels, numnodes)
% Function that processes the output of the Shen method into the matrix

% Finds the number of communities, and makes a matrix of zeros.
Output = zeros(numnodes, max([NodeLabels{:}]));

for i = 1:numnodes
    % equally separates the nodes into communities
    Output(i, NodeLabels{i}) = 1/length(NodeLabels{i});
end