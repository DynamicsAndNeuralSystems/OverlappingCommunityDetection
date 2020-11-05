function [Output] = call_Jerry(Mat, numnodes, numIters, Threshold, IsBlind, ListeningRule)
% Jerry
% Input is Adjacency matrix
% Check for inputs
if nargin < 3 || isempty(numIters)
    numIters = 120;
end
if nargin < 4 || isempty(Threshold)
    Threshold = 0.09;
end
if nargin < 5 || isempty(IsBlind)
    IsBlind = 1;
end
if nargin < 6 || isempty(ListeningRule)
%     ListeningRule = 'majority';
    ListeningRule = 'probabilistic';
end

% Runs the Jerry method
[Jerry] = run_Jerry(Mat, numIters, Threshold, IsBlind, ListeningRule, numnodes);

% Processes the output
[Jerry_final] = process_Jerry(Jerry, numnodes); 

% Places the structure data in
Output = struct('Name', 'Jerry', 'Result', Jerry_final);
