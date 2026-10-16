function M = communitiesToMatrix(communityNodes,numNodes)
% Converts a cell array of community node lists (one cell per community) into
% a numNodes x numCommunities membership matrix, tidied by tidyCommunities
% (overlapping nodes split their membership equally across communities).
%-------------------------------------------------------------------------------

numComms = length(communityNodes);
M = zeros(numNodes,numComms);
for c = 1:numComms
    M(communityNodes{c},c) = 1;
end
M = tidyCommunities(M);

end
