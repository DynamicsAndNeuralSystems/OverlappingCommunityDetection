function M = process_Gopalan(outFile, numNodes)
% Processes the text-file output (communities.txt) of the Gopalan method:
% one community per line, listing the (original) ids of its nodes.

communityNodes = readCommunityLines(outFile);
M = communitiesToMatrix(communityNodes, numNodes);

end
