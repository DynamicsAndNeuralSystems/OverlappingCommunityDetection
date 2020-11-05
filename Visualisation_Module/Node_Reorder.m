function Output = Node_Reorder(CommMat)
% Reorders the matrix so that 1 is on top.

numComms = size(CommMat, 2); % Finds the number of communities
Output = CommMat; % Keeps the matrix for the output

% CommMat = ~(CommMat == 0); % Creates a logical of every connection

Sizes = sum(CommMat, 1); % Finds the sizes of each community

[~, I] = sort(Sizes, 'descend'); % Sorts them in terms of size

Output = Output(:, I); % Reorders communities

end
