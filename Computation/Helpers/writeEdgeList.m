function writeEdgeList(fileName,edgeList,includeWeights)
% Writes an edge list [node1, node2, weight] to a tab-separated text file.
%-------------------------------------------------------------------------------

if nargin < 3
    includeWeights = true;
end
fid = fopen(fileName,'w');
if fid < 0
    error('Could not open %s for writing',fileName);
end
if includeWeights
    fprintf(fid,'%u\t%u\t%.10g\n',edgeList(:,1:3)');
else
    fprintf(fid,'%u\t%u\n',edgeList(:,1:2)');
end
fclose(fid);

end
