function outFile = run_Link(Undir, sourceCodePath, workDir)
% Runs the link communities command (weighted, tab-delimited input) inside the
% existing, writable directory workDir. Returns the path to the comm2nodes file.
% The python executable is taken from the environment variable OCDA_PYTHON if
% set, otherwise 'python3'.

pythonExe = getenv('OCDA_PYTHON');
if isempty(pythonExe)
    pythonExe = 'python3';
end
scriptFile = fullfile(sourceCodePath, 'Link', 'link_clustering.py');
if ~isfile(scriptFile)
    error('Cannot find %s', scriptFile);
end

netFile = fullfile(workDir, 'Link.txt');
writeEdgeList(netFile, Undir, true);

% -w: use weights (writeEdgeList is tab-delimited, the script's default)
runSystemCommand(sprintf('"%s" "%s" -w "%s"', pythonExe, scriptFile, netFile));

% Output is named from the input basename (Link_maxS<S>_maxD<D>.comm2nodes.txt)
outFiles = dir(fullfile(workDir, '*comm2nodes.txt'));
if numel(outFiles) ~= 1
    error('Expected one comm2nodes output file in %s, found %d', workDir, numel(outFiles));
end
outFile = fullfile(workDir, outFiles(1).name);

end
