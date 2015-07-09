function [Output] = process_Clauset(Input, numnodes)
% Function that processes the data output of Clauset into a nice matrix for
% the visualiser to use

% Preallocates the matrix
Output = zeros(numnodes, 1);

% Places 1s where the nodes belong in the community
for i = 1:numnodes
    Output(i, Input(i)) = 1;
end