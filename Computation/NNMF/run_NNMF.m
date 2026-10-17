function NNMF = run_NNMF(AdjMat,thresh)
% Runs the NNMF algorithm (commDetNMF), then thresholds the membership matrix
% independently at each threshold. Returns a cell with one matrix per threshold.
%-------------------------------------------------------------------------------

P = commDetNMF(AdjMat);

NNMF = cell(length(thresh),1);
for i = 1:length(thresh)
    % Remove memberships below threshold; unassigned nodes get zero rows
    NNMF{i} = tidyCommunities(P.*(P > thresh(i)));
end

end
