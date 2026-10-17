function I = Node_Sorter(benchComms)
% Sorts nodes according to community.
% Nodes are grouped by the first benchmark community they belong to; nodes
% belonging to no community (all-zero or NaN rows) are appended at the end.
% The output is a full permutation (column vector) of 1:numNodes.
%-------------------------------------------------------------------------------

numNodes = size(benchComms, 1);

benchComms = ~(benchComms == 0 | isnan(benchComms)); % Creates a logical of every connection

% Preallocates the matrix for the sorting
I = zeros(numNodes, 1);

% Counter for putting things in
counter = 0;

for i = 1:size(benchComms, 2)
    nodesHere = find(benchComms(:, i));
    % Places values into the sort
    I(counter+1:counter+length(nodesHere)) = nodesHere;

    % Changes the counter
    counter = counter + length(nodesHere);

    % Turns those lines that have been used to 0.
    benchComms(nodesHere, :) = false;
end

% Append nodes that are in no community:
used = false(numNodes, 1);
used(I(1:counter)) = true;
I(counter+1:end) = find(~used);

end
