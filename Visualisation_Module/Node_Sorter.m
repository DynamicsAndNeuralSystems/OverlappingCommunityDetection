function [I] = Node_Sorter(BenchComms)
% Function that sorts nodes according to community.

% Preallocates the matrix for the sorting
I = zeros(size(BenchComms, 1), 1);

BenchComms = ~(BenchComms == 0); % Creates a logical of every connection

% Counter for putting things in
counter = 0;

for i = 1:size(BenchComms, 2)
    % Places values into the sort
    I(counter+1:counter+sum(BenchComms(:, i))) = find(BenchComms(:, i) == 1);
    
    % Changes the counter
    counter = counter + sum(BenchComms(:, i));
    
    % Turns those lines that have been used to 0.
    BenchComms(BenchComms(:, i) == 1, :) = 0;
    
end
