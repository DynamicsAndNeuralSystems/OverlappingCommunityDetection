function DirList = Undir2Direct(Undir)
% Converts an undirected edge list to a directed list (each edge in both
% directions; self-loops are not duplicated), sorted by from-node then to-node.
%-------------------------------------------------------------------------------

isSelf = (Undir(:,1) == Undir(:,2));
DirList = [Undir; Undir(~isSelf,[2 1 3])];
DirList = sortrows(DirList,[1 2]);

end
