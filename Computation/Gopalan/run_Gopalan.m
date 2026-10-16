function outFile = run_Gopalan(Undir, numNodes, maxIter, sourceCodePath, workDir)
% Runs the Gopalan method (SVINET) on the network, working inside workDir
% (an existing, writable directory). Returns the path to the communities file
% of the final (link-sampling) run.

svinet = fullfile(sourceCodePath, 'Gopalan', 'src', 'bin', 'svinet');
if ~isfile(svinet)
    error(['Cannot find the SVINET executable:\n  %s\nBuild it with ' ...
        '"./configure && make" in %s (requires the GNU Scientific Library, GSL).'], ...
        svinet, fullfile(sourceCodePath, 'Gopalan'));
end
[~, attr] = fileattrib(svinet);
if ~attr.UserExecute
    error('SVINET (%s) is not executable: run chmod +x on it, or rebuild it', svinet);
end

netFile = fullfile(workDir, 'network.txt');
writeEdgeList(netFile, Undir, false);

% svinet writes its output directories into the current directory, so work in workDir
origDir = pwd;
restoreDir = onCleanup(@()cd(origDir)); %#ok<NASGU>
cd(workDir);

% Find the approximate number of communities (output: n<N>-k<K>-<label>-findk/)
runSystemCommand(sprintf('"%s" -file "%s" -n %d -k %d -label temp -findk', ...
                                svinet, netFile, numNodes, numNodes));
findkFile = fullfile(workDir, sprintf('n%d-k%d-temp-findk', numNodes, numNodes), ...
                                'communities.txt');
numComms = length(readCommunityLines(findkFile));
if numComms == 0
    error('SVINET found no communities');
end

% Run the algorithm (output: n<N>-k<K>-<label>-linksampling/)
runSystemCommand(sprintf(['"%s" -file "%s" -n %d -k %d -label Gopalan ' ...
        '-link-sampling -max-iterations %d'], svinet, netFile, numNodes, numComms, maxIter));
outFile = fullfile(workDir, sprintf('n%d-k%d-Gopalan-linksampling', numNodes, numComms), ...
                                'communities.txt');
if ~isfile(outFile)
    error('SVINET did not produce the expected output file %s', outFile);
end

end
