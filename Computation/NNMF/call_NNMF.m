function Output = call_NNMF(Mat, thresh)
% NNMF
% 
%---INPUTS:
% Mat is the adjacency matrix
% thresh are the threshold(s) for denoting membership to a community
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check Inputs
%-------------------------------------------------------------------------------
% Thresholds (can be multiple)
if nargin < 2 || isempty(thresh)
    thresh = [0,0.01,0.1];
end
%-------------------------------------------------------------------------------

NNMF_final = run_NNMF(Mat, thresh); % Runs the NNMF method

NNMF_final = NNMF_final{1}; % Gets rid of the cells

% Gets rid of all communities that no longer exist due to thresholds
NNMF_final = NNMF_final(:, sum(NNMF_final, 1) ~= 0);

% No need for futher processing: the output is the required matrix format

% Puts the structural data in
Output = struct('Name', 'NNMF', 'Threshold', thresh, ...
                    'Result', NNMF_final);

end
