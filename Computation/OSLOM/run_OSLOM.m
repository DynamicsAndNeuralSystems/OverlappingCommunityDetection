function [] = run_OSLOM(Undir, numIters, Tol)
% Run the OSLOM algorithm (in Linux commandline) from Matlab
%-------------------------------------------------------------------------------

filePath = './Computation_Module/SourceCode/OSLOM'; % File path from module to OSLOM code
fileBack = '../../../'; % File path back

cd(filePath); % Goes into the pathway of the code

fid = fopen('tempData.txt', 'w'); % Creates a text file
fprintf(fid, '%g\t%g\t%f\n', Undir'); % Saves the sparse matrix as a textfile
fclose(fid);

for x = Tol
    % Runs the command to start the OSLOM algorithm
    system(sprintf('./oslom_undir -f tempData.txt -w -r %g -t %f', numIters, x));
    % Moves and renames the important data to a nicer location
    system(sprintf('mv tempData.txt_oslo_files/tp %s/OSLOM_tol_%g.txt', fileBack, x));
    % Removes unneeded files
    system('rm -r tempData.txt_oslo_files');
    
end

system('rm tempData.txt'); % Deletes the text file
cd(fileBack); % Moves back to the first directory

end
