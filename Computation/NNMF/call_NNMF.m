function Output = call_NNMF(Mat,thresh)
% NNMF (nonnegative matrix factorization)
%
%---INPUTS:
% Mat: adjacency matrix
% thresh: scalar threshold for membership of a community (default 0.1)
%-------------------------------------------------------------------------------

if nargin < 2 || isempty(thresh)
    thresh = 0.1;
end

NNMF_final = run_NNMF(Mat,thresh);

Output = struct('Name','NNMF','Threshold',thresh,'Result',NNMF_final{1});

end
