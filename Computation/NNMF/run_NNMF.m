function NNMF = run_NNMF(AdjMat, thresh)
% Runs the NNMF algorithm.
% Can be adapted to include a threshold
%-------------------------------------------------------------------------------

NNMF = cell(length(thresh), 1);

% Performs NNMF:
temp = commDetNMF(AdjMat);

for i = thresh
    temp = temp.*(temp > i); % Deletes values that are below the threshold    
    for x = 1:length(AdjMat)
        % Extends all other values to make up for the missing values
        NNMF{thresh==i}(x, :) = temp(x, :)/sum(temp(x, :));
    end
end

end
