function [] = run_Link(Undir)
% Function that runs the link communities command

temp = Undir(:, 1:2); % Creates the matrix for the input

% Goes into the directory required
cd Modules/Computation_Module/Computation/Link;

fid = fopen('Link.txt', 'w'); % Creates text file
fprintf(fid, '%g\t%g\n', temp'); % Places data within file
fclose(fid);

