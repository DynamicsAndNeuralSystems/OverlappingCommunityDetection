function Cover = ConvertNodeLabelsToCover(NodeLabels)

try
    NumComms = length(unique(vertcat(NodeLabels{:})));
catch
    NumComms = length(unique(horzcat(NodeLabels{:})));
end
NumNodes = length(NodeLabels);

Cover = zeros(NumNodes,NumComms);

for i = 1:NumComms
    Cover(:,i) = cellfun(@(x)ismember(i,x),NodeLabels);
end

end