function I = Node_Sorter(benchComms)
% Sorts nodes according to community.
%-------------------------------------------------------------------------------

% Preallocates the matrix for the sorting
I = zeros(size(benchComms, 1), 1);

benchComms = ~(benchComms == 0); % Creates a logical of every connection

% Counter for putting things in
counter = 0;

for i = 1:size(benchComms, 2)
    % Places values into the sort
    I(counter+1:counter+sum(benchComms(:, i))) = find(benchComms(:, i) == 1);
    
    % Changes the counter
    counter = counter + sum(benchComms(:, i));
    
    % Turns those lines that have been used to 0.
    benchComms(benchComms(:, i) == 1, :) = 0;
    
end

end
