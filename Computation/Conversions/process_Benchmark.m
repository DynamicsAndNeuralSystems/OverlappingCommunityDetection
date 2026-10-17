function Output = process_Benchmark(benchFilename,numNodes)
% Processes a benchmark community file into a numNodes x numComms membership
% matrix. Each (non-comment) line is "node comm1 comm2 ...": the first value is
% the node id and the rest are the communities it belongs to (nodes need not be
% listed in order, and trailing whitespace is irrelevant). Nodes in several
% communities split their membership equally.
%-------------------------------------------------------------------------------

lineValues = readCommunityLines(benchFilename);

numComms = 0;
for i = 1:length(lineValues)
    numComms = max([numComms;lineValues{i}(2:end)]);
end

Output = zeros(numNodes,numComms);
for i = 1:length(lineValues)
    v = lineValues{i};
    Output(v(1),v(2:end)) = 1;
end
Output = tidyCommunities(Output);

end
