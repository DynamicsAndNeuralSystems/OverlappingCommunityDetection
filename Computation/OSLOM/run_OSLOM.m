function tpFile = run_OSLOM(Undir,numIters,tol,sourceCodePath,workDir)
% Runs the OSLOM binary (undirected, weighted) on an edge list at a single
% tolerance, working entirely inside the temporary directory workDir.
% Returns the path of the output file listing the modules found.
%-------------------------------------------------------------------------------

oslomBinary = fullfile(sourceCodePath,'OSLOM','oslom_undir');
if ~isfile(oslomBinary)
    error('OSLOM binary not found at %s: run SourceCode/OSLOM/compile_all.sh',oslomBinary);
end

% Write the edge list:
inputName = 'tempData.txt';
writeEdgeList(fullfile(workDir,inputName),Undir,true);

% OSLOM overflows a fixed-size path buffer on long file paths (as in macOS
% temp directories), so run from within workDir using a short relative name.
% Outputs go to <inputName>_oslo_files, next to the input.
originalDir = pwd;
returnToDir = onCleanup(@()cd(originalDir)); %#ok<NASGU>
cd(workDir);
fprintf(1,'Running OSLOM at tolerance %g\n',tol);
runSystemCommand(sprintf('"%s" -f %s -w -r %g -t %f',oslomBinary,inputName,numIters,tol));

tpFile = fullfile(workDir,[inputName,'_oslo_files'],'tp');

end
