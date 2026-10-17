function Output = call_infomap(DirList,numNodes,sourceCodePath)
% Infomap (overlapping modules)
%
%---INPUTS:
% DirList: directed edge list [from, to, weight]
% numNodes: number of nodes
% sourceCodePath: path to SourceCode (containing Infomap/Infomap)
%-------------------------------------------------------------------------------

if nargin < 3 || isempty(sourceCodePath)
    sourceCodePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'SourceCode');
end

[workDir,cleanupObj] = makeTempDir(); %#ok<ASGLU>
cluFile = run_infomap(DirList,sourceCodePath,workDir);
Infomap_final = process_infomap(cluFile,numNodes);

Output = struct('Name','Infomap','Result',Infomap_final);

end
