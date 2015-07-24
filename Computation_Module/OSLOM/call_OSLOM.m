function [Output] = call_OSLOM(Undir, numnodes, numIter, Tol)
%% OSLOM
% Input is a sparse matrix (undirected)
% Check inputs
if nargin < 3 || isempty(numIter)
    numIter = 100; % Number of iterations within the algorithm
end
if nargin < 4 || isempty(Tol)
    Tol = 0.5; % Range of tolerances
end

% Runs the algorithm for OSLOM
run_OSLOM(Undir, numIter, Tol);

% Processes the textfile outputs
[OSLOM_final] = process_OSLOM(Tol, numnodes);

    % Puts the structural data in
    Output = struct('Name', 'OSLOM', 'Threshold', Tol, ...
        'Result', OSLOM_final{1});