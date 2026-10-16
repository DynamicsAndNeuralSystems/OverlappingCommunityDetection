function clu = run_infomap(DirList,sourceCodePath,workDir)
% Runs the Infomap binary (with overlapping modules) on a directed edge list,
% working entirely inside the temporary directory workDir.
% Returns the path of the output .clu file.
%-------------------------------------------------------------------------------

infomapBinary = fullfile(sourceCodePath,'Infomap','Infomap');
if ~isfile(infomapBinary)
    error('Infomap binary not found at %s: build it with make in SourceCode/Infomap',infomapBinary);
end

inputFile = fullfile(workDir,'networki.txt');
writeEdgeList(inputFile,DirList,true);
outDir = fullfile(workDir,'output');
mkdir(outDir);

runSystemCommand(sprintf('"%s" "%s" "%s" --overlapping --clu',infomapBinary,inputFile,outDir));
clu = fullfile(outDir,'networki.clu');

end
