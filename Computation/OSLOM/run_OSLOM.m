function run_OSLOM(Undir, numIters, Tol, sourceCodePath, outputPath)
% Run the OSLOM algorithm (in Linux commandline) from Matlab
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Write data into the OSLOM directory
%-------------------------------------------------------------------------------
% pathBase = GiveMeBaseDir();
% pathBase = fileparts(fileparts(pwd));
% filePath = fullfile('SourceCode','OSLOM'); % File path from module to OSLOM code
cd(fullfile(sourceCodePath, 'OSLOM')); % Goes into the pathway of the code

%-------------------------------------------------------------------------------
% Write the sparse matrix out to a text file
%-------------------------------------------------------------------------------
fid = fopen('tempData.txt', 'w');
fprintf(fid,'%g\t%g\t%f\n', Undir');
fclose(fid);

%-------------------------------------------------------------------------------
% Run OSLOM
%-------------------------------------------------------------------------------
numTol = length(Tol);
fprintf(1,'Running OSLOM on the network across %u tolerances\n',numTol);
for t = 1:numTol
    tol = Tol(t);
    % Run the command to start the OSLOM algorithm:
    system(sprintf('./oslom_undir -f tempData.txt -w -r %g -t %f', numIters, tol));
    % Move/rename the important data to base location:
    system(sprintf('mv tempData.txt_oslo_files/tp %s/OSLOM_tol_%g.txt', outputPath, tol));
    % Remove unneeded files:
    system('rm -r tempData.txt_oslo_files');
end

%-------------------------------------------------------------------------------
% Clean up
%-------------------------------------------------------------------------------
% Delete the temporary data file:
system('rm tempData.txt'); 
% Return to base directory:
% cd(pathBase);

end
