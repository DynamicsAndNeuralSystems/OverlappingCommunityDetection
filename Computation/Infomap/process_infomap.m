function Output = process_infomap(cluFile,numNodes)
% Processes an Infomap .clu file (lines "node module flow", '#' comments) into
% a numNodes x numModules membership matrix. With --overlapping, a node may
% appear in several modules (its membership is split equally).
%-------------------------------------------------------------------------------

lines = readCommunityLines(cluFile);
A = cell2mat(cellfun(@(x)x(1:2)',lines,'UniformOutput',false)); % [node, module]
numComms = max(A(:,2));
Output = zeros(numNodes,numComms);
Output(sub2ind([numNodes,numComms],A(:,1),A(:,2))) = 1;
Output = tidyCommunities(Output);

end
