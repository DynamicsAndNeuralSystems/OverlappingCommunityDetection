function [Output] = process_Clauset(Input,numNodes)
% Function that processes the data output of Clauset into a nice matrix for
% the visualizer to use (hard partition: one block per node; empty blocks dropped)

% Preallocates the matrix
Output = zeros(numNodes,max(Input));

% Places 1s where the nodes belong in the community
for i = 1:numNodes
    Output(i,Input(i)) = 1;
end

Output = tidyCommunities(Output);

end
