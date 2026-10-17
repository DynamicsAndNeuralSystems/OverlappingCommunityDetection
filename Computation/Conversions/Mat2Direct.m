function DirList = Mat2Direct(Mat,numNodes)
% Converts an adjacency matrix to a directed edge list [from, to, weight],
% ordered by from-node and then to-node.
%-------------------------------------------------------------------------------

if nargin < 2
    numNodes = size(Mat,1);
end
Mat = Mat(1:numNodes,1:numNodes);

% find on the transpose returns entries in row-major order of Mat:
[to,from] = find(Mat');
DirList = [from,to,Mat(sub2ind(size(Mat),from,to))];

end
