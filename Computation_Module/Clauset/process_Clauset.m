function [Output] = process_Clauset(Input,numNodes)
% Function that processes the data output of Clauset into a nice matrix for
% the visualizer to use

% Preallocates the matrix
Output = zeros(numNodes, 1);

% Places 1s where the nodes belong in the community
for i = 1:numNodes
    Output(i, Input(i)) = 1;
end

end
