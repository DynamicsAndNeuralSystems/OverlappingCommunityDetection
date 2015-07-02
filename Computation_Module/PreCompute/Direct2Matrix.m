function [Mat] = Direct2Matrix(DirList, numnodes)
% Converts the directed data into a matrix
Mat = zeros(numnodes, numnodes); % Allocates a 0 matrix for creation

for x = 1:size(DirList, 1)
    Mat(DirList(x, 1), DirList(x, 2)) = DirList(x, 3); % Sets each value to the power
end