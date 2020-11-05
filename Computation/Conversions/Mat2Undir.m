function [Undir] = Mat2Undir(Mat)
% Function that converts an adjacency matrix to an undirected list

% Finds the indices of every value that is non-zero in both sides of the
% matrix (both below and above the self connection edges)
MatUnd = ((Mat~=0) + (Mat~=Mat')) ~= 0;

Undir= [];

% For all nodes
for x = 1:size(Mat, 2)
    % For all nodes above (x,y)
    for y = x:size(Mat,2)
        % If there is a value
        if MatUnd(x, y) == 1
            % Put this edge into the undirected list
            Undir(end+1, :) = [x, y, (Mat(x, y) + Mat(y, x))/2];
        end
    end
end