function [NumCommMat, NumOverlap, NumHomeles] = Plot_Plan(NumCommMat, NumOverlap, NumHomeles, NodeLabel)
% Function that finds then adds the number of communities, number of
% overlapping nodes and number of homeless nodes to the original matrices,
% to eventually show it on a plot.

% Creation of the Number of communities matrix
NumCommMat(1, end + 1) = max([NodeLabel{:}]);

% Creation of the overlapping node communities
NumOverlap(1, end + 1) = sum(cellfun('length', NodeLabel) > 1);

% Creation of the homeless nodes matrix
NumHomeles(1, end + 1) = sum(cellfun('length', NodeLabel) < 1);