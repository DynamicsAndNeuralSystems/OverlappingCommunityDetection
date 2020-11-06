function Cover = ConvertNodeLabelsToCover(nodeLabels)

try
    numComms = length(unique(vertcat(nodeLabels{:})));
catch
    numComms = length(unique(horzcat(nodeLabels{:})));
end
numNodes = length(nodeLabels);

Cover = zeros(numNodes,numComms);

for i = 1:numComms
    Cover(:,i) = cellfun(@(x)ismember(i,x),nodeLabels);
end

end
