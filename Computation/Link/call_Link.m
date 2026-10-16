function Output = call_Link(Undir, numNodes, prec, sourceCodePath)
% Link communities (Ahn, Bagrow & Lehmann, 2010)
% Input is undirected list [node1, node2, weight]; node ids in 1..numNodes
% sourceCodePath is the absolute path to the toolbox's SourceCode directory.

if nargin < 3 || isempty(prec)
    prec = 0.01; % Proportion of lowest-weight links to cut before computation
end
if nargin < 4 || isempty(sourceCodePath)
    error('Specify the path to the SourceCode directory');
end

% Gets rid of the weakest links (by the proportion prec)
Input = removeWeakestLinks(Undir, prec);

% Runs the link communities method in a temporary directory
[workDir, cleanupObj] = makeTempDir(); %#ok<ASGLU>
outFile = run_Link(Input, sourceCodePath, workDir);

% Processes the text-file output
Link_final = process_Link(outFile, numNodes);

% Puts the structural data in
Output = struct('Name', 'Link', 'Percent_removed', prec*100, 'Result', Link_final);

end
