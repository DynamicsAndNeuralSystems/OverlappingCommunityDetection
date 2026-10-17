function Mat = Direct2Matrix(DirList,numNodes)
% Converts a directed edge list [from, to, weight] into an adjacency matrix.
% If an edge is repeated, the last occurrence sets the weight.
%-------------------------------------------------------------------------------

Mat = zeros(numNodes,numNodes);
Mat(sub2ind([numNodes,numNodes],DirList(:,1),DirList(:,2))) = DirList(:,3);

end
