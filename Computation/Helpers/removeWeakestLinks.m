function Undir = removeWeakestLinks(Undir,prop)
% Removes the weakest proportion, prop (in [0,1)), of links from an edge list
% [node1, node2, weight]. Links are ranked by weight (ties broken by list
% order), so exactly round(prop*numLinks) links are removed, including on
% unweighted networks where all weights tie.
%-------------------------------------------------------------------------------

if prop <= 0
    return
end
numLinks = size(Undir,1);
[~,ix] = sort(Undir(:,3),'ascend');
Undir(ix(1:round(prop*numLinks)),:) = [];

end
