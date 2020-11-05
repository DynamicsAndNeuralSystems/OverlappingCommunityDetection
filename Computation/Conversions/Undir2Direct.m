function [DirList] = Undir2Direct(Undir)
% Function that converts an undirected list to a directed list

DirList = [Undir; Undir(:, [2 1 3])]; % Doubles the list

DirList = sortrows(DirList, [1 2]); % Sorts everything in the first row, then in the second