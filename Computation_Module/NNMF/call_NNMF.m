function [Output] = call_NNMF(Mat, thresh)
% NNMF
% Input is Adjacency matrix
% Check inputs
if nargin < 2 || isempty(thresh)
    thresh = [0 0.01 0.1]; % Thresholds, can be multiple thresholds.
end

NNMF_final = run_NNMF(Mat, thresh); % Runs the NNMF method

% This has no need for processing, as the output is the required matrix

    % Puts the structural data in
    Output = struct('Name', 'NNMF', 'Threshold', thresh, ...
        'Result', NNMF_final{1});
