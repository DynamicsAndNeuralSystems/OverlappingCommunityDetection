function Output = call_OSLOM(Undir,numNodes,numIter,tol,sourceCodePath)
% OSLOM (Order Statistics Local Optimization Method)
%
%---INPUTS:
% Undir: undirected edge list [node1, node2, weight]
% numNodes: number of nodes
% numIter: number of iterations of the algorithm (default 100)
% tol: significance tolerance (scalar, default 0.5)
% sourceCodePath: path to SourceCode (containing OSLOM/oslom_undir)
%-------------------------------------------------------------------------------

if nargin < 3 || isempty(numIter)
    numIter = 100;
end
if nargin < 4 || isempty(tol)
    tol = 0.5;
end
if nargin < 5 || isempty(sourceCodePath)
    sourceCodePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'SourceCode');
end

% Run in a temporary directory (deleted on return):
[workDir,cleanupObj] = makeTempDir(); %#ok<ASGLU>
tpFile = run_OSLOM(Undir,numIter,tol,sourceCodePath,workDir);
OSLOM_final = process_OSLOM(tpFile,numNodes);

Output = struct('Name','OSLOM','Threshold',tol,'Result',OSLOM_final);

end
