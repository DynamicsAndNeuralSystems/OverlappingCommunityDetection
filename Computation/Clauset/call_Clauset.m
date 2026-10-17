function Output = call_Clauset(DirList, numNodes, kRange)
% Clauset - weighted stochastic block model (WSBM; Aicher, Jacobs & Clauset)
% Uses a weighted directed list of connections as an input.
% NOTE: WSBM gives hard (non-overlapping) partitions: each node is in exactly one block.
% The number of blocks k is chosen from kRange (default 2:8) as the value maximizing
% the log-evidence of the fitted model (Model.Para.LogEvidence).
%---Output fields: Name, k (chosen), Result (node-by-community matrix),
% kRange and LogEvidence (per-k log-evidence).
%-------------------------------------------------------------------------------

if nargin < 3 || isempty(kRange)
    kRange = 2:8;
end

numK = length(kRange);
logEvidence = nan(numK,1);
allLabels = cell(numK,1);
for i = 1:numK
    [allLabels{i},Model] = run_Clauset(DirList,kRange(i));
    logEvidence(i) = Model.Para.LogEvidence;
    fprintf(1,'WSBM with k = %u: log-evidence = %g\n',kRange(i),logEvidence(i));
end

[~,bestInd] = max(logEvidence);
chosenK = kRange(bestInd);
fprintf(1,'Chose k = %u (maximum log-evidence)\n',chosenK);

Clauset_final = process_Clauset(allLabels{bestInd},numNodes); % Processes the output

% Putting it in the structure
Output = struct('Name','Clauset','k',chosenK,'Result',Clauset_final, ...
                'kRange',kRange,'LogEvidence',logEvidence);

end
