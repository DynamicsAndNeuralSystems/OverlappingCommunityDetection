function M = process_Link(outFile, numNodes)
% Processes the comm2nodes text-file output of the link communities method:
% one community per line, the first value is the community id, the rest are node ids.

lineValues = readCommunityLines(outFile);
communityNodes = cellfun(@(x)x(2:end), lineValues, 'UniformOutput', false);
M = communitiesToMatrix(communityNodes, numNodes);

end
