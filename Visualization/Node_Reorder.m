function Output = Node_Reorder(CommMat)
% Reorders the matrix so that the largest community is first (on top).
% NaN entries are treated as 0.

CommMat(isnan(CommMat)) = 0;

Sizes = sum(CommMat, 1); % Finds the sizes of each community

[~, I] = sort(Sizes, 'descend'); % Sorts them in terms of size

Output = CommMat(:, I); % Reorders communities

end
