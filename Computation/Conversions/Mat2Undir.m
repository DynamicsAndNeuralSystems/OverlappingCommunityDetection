function Undir = Mat2Undir(Mat)
% Converts an adjacency matrix to an undirected edge list [node1, node2, weight]
% with node1 <= node2 (self-loops included), ordered by node1 then node2.
% Weights are symmetrized as (Aij + Aji)/2, so a one-directional edge
% (only Aij nonzero) receives half its weight.
%-------------------------------------------------------------------------------

% Edges present in either direction, upper triangle (including diagonal):
isEdge = triu((Mat ~= 0) | (Mat' ~= 0));
[col,row] = find(isEdge'); % row-major order
Undir = [row,col,(Mat(sub2ind(size(Mat),row,col)) + Mat(sub2ind(size(Mat),col,row)))/2];

end
