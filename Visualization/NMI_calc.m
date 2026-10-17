function NMI = NMI_calc(Input1, Input2)
% Calculates the normalized similarity of two community assignments.
% **[NON-OVERLAPPING]**
% Inputs are numNodes x numComms membership matrices describing hard
% partitions (at most one nonzero entry per row). Overlapping inputs give an
% error (use ENMI_calc for overlapping covers).
% Nodes unassigned (all-zero or NaN row) in either input are excluded from the
% calculation (with a warning).
% For hard partitions the result equals NMI_calc_vecs on the community labels.
%-------------------------------------------------------------------------------

% Binarize (NaN counts as unassigned):
Comm1 = (Input1 ~= 0) & ~isnan(Input1);
Comm2 = (Input2 ~= 0) & ~isnan(Input2);

if size(Comm1,1)~=size(Comm2,1)
    error('Inputs must have the same number of nodes (rows)');
end
if any(sum(Comm1,2) > 1) || any(sum(Comm2,2) > 1)
    error('NMI_calc is only valid for non-overlapping (hard) partitions: some node has more than one community. Use ENMI_calc for overlapping covers.');
end

% Exclude unassigned nodes
assigned = (sum(Comm1,2)==1) & (sum(Comm2,2)==1);
if ~all(assigned)
    warning('%u unassigned node(s) excluded from the NMI calculation',sum(~assigned));
end
if ~any(assigned)
    NMI = 0;
    return
end

% Community label of each assigned node (first, and only, nonzero column)
[~, labels1] = max(Comm1(assigned,:), [], 2);
[~, labels2] = max(Comm2(assigned,:), [], 2);

NMI = NMI_calc_vecs(labels1, labels2);

end
