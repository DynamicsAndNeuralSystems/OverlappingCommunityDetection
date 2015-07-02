function[Matrix] = NodeLabel2Mat(NodeLabel)
% Function that takes a cell as an input and outputs a matrix used for
% rectangles in the visualiser
% Brandon Lam, 23-07-2014

% Preallocates
Matrix = zeros(size(NodeLabel, 1), max([NodeLabel{:}]));

% Loop converts it to a matrix
for i = 1:size(NodeLabel, 1)
    % Adds 1/[no. of comms] as a value for every place
    Matrix(i, NodeLabel{i}) = 1/size(NodeLabel{i}, 2);
end