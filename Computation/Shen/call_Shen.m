function Output = call_Shen(Mat, numNodes, cliqueSize, cwThresh, lThresh)
% OCDA method of Shen et al. (2009):
% 
% Shen et al. (2009). Quantifying and identifying the overlapping community
% structure in networks. Journal of Statistical Mechanics: Theory and Experiment.
% http://stacks.iop.org/1742-5468/2009/i=07/a=P07042
% 
%---Inputs:
% Adjacency matrix, Mat
% Number of nodes, numNodes
% Clique size, k
% Clique weight threshold, cwThresh
% link threshold, lThresh
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs
%-------------------------------------------------------------------------------
if nargin < 3 || isempty(cliqueSize)
    cliqueSize = 7;
end
if nargin < 4 || isempty(cwThresh)
    cwThresh = 1.3; % Clique weight threshold
end
if nargin < 5 || isempty(lThresh)
    lThresh = 0; % Threshold on links, 0 because we want to preserve data
end
%-------------------------------------------------------------------------------

% Runs the function for the Shen method
Shen = run_Shen(Mat, cliqueSize, cwThresh, lThresh);

% Processes the data into the matrix
Shen_final = process_Shen(Shen, numNodes);

% Places the structure data in
Output = struct('Name','Shen','Result',Shen_final);

end
