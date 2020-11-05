function [] = run_Link(Undir)
% Function that runs the link communities command

% Goes into the directory required
cd Computation_Module/SourceCode/Link;

fid = fopen('Link.txt', 'w'); % Creates text file
fprintf(fid, '%g\t%g\t%f\n', Undir'); % Places data within file
fclose(fid);

system('./link_clustering.py Link.txt');
system('mv *comm2nodes.txt ../../../../Link.txt');
system('rm *.txt');

cd ../../../;
