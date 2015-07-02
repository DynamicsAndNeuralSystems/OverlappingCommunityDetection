function [Undir] = Direct2Undir(DirList)
% Function that converts a directed list to an undirected list

% Cuts out all connections that are repeated (if the first node is lower
% in number than the second in the connection)
Undir = DirList(DirList(:, 1)<DirList(:, 2), :);