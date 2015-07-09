function [Output] = call_Shen(Mat, numnodes, cSize, cwThresh, lThresh)
% Shen
% Input is Adjacency matrix
% Check inputs
if nargin < 3 || isempty(cSize)
    cSize = 7; % Clique size
end
if nargin < 4 || isempty(cwThresh)
    cwThresh = 1.3; % Clique weight threshold
end
if nargin < 5 || isempty(lThresh)
    lThresh = 0; % Threshold on links, 0 because we want to preserve data
end

% Runs the function for the Shen method
[Shen] = run_Shen(Mat, cSize, cwThresh, lThresh);

% Processes the data into the matrix
[Shen_final] = process_Shen(Shen, numnodes);

% Places the structure data in
Output = struct('Name', 'Shen', 'Result', Shen_final);