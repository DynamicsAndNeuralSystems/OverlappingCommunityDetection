function[sorted_comm, Index] = comm_sort(Comm_Cell)
% Sorts out in terms of community, and outputs the index
% Brandon Lam, 11-07-2014

% Preallocates a cell matrix
sorted_comm = cell(1, length(Comm_Cell)); 

for i = 1:length(Comm_Cell) % For all nodes
    % Strip the communities off until there is only one for each node
    sorted_comm{i} = Comm_Cell{i}(1); 
    % This will be used for sorting
end

% The actual sort - this is in order of communities
[sorted_comm, Index] = sort(cell2mat(sorted_comm)); 

end