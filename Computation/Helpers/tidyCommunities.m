function M = tidyCommunities(M)
% Tidies a node-by-community membership matrix into the toolbox's standard format:
% - NaN entries are set to 0
% - each row with any membership is normalized to sum to 1
%   (rows of unassigned nodes stay all-zero rather than becoming 0/0 = NaN)
% - empty communities (all-zero columns) are removed
%-------------------------------------------------------------------------------

M(isnan(M)) = 0;
rowSums = sum(M,2);
isAssigned = rowSums > 0;
M(isAssigned,:) = M(isAssigned,:)./rowSums(isAssigned);
M = M(:,any(M > 0,1));

end
