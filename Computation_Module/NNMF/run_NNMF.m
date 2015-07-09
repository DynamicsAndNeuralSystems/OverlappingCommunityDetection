function [NNMF] = run_NNMF(AdjMat, thresh)
% Function that calls the NNMF algorithm.
% Can be adapted to include a threshold
NNMF = cell(length(thresh), 1);

temp = commDetNMF(AdjMat); % Does the NNMF calculation

for i = thresh
    temp = temp.*(temp > i); % Deletes values that are below the threshold
    
    for x = 1:length(AdjMat)
        % Extends all other values to make up for the missing values
        NNMF{thresh==i}(x, :) = temp(x, :)/sum(temp(x, :));
    end
end