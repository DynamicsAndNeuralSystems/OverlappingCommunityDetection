function Output = call_Gopalan(Undir, numNodes, maxIter, prec, sourceCodePath)
% Gopalan (SVINET)
% Input is either undirected or directed list, but it will treat them as
% undirected. Node ids must be integers in 1..numNodes.
% sourceCodePath is the absolute path to the toolbox's SourceCode directory.

%-------------------------------------------------------------------------------
%% Check for inputs
%-------------------------------------------------------------------------------
if nargin < 3 || isempty(maxIter)
    maxIter = 1000; % Maximum iterations done within the algorithm
end
if nargin < 4 || isempty(prec)
    prec = 0.1; % Proportion of lowest-weight links to cut before computation
end
if nargin < 5 || isempty(sourceCodePath)
    error('Specify the path to the SourceCode directory');
end

% Gets rid of the weakest links (by the proportion prec)
Input = removeWeakestLinks(Undir, prec);

% Runs the Gopalan method in a temporary directory
[workDir, cleanupObj] = makeTempDir(); %#ok<ASGLU>
outFile = run_Gopalan(Input, numNodes, maxIter, sourceCodePath, workDir);

% Processes the output (unassigned nodes get zero rows; empty communities removed)
Gopalan_final = process_Gopalan(outFile, numNodes);

% Puts the structural data in
Output = struct('Name', 'Gopalan', 'Percent_removed', prec*100, ...
    'Max_Iterations', maxIter, 'Result', Gopalan_final);

end
