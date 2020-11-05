function [] = run_Gopalan(Undir, numnodes, maxIter)
% Function that runs the Gopalan method for the network

temp = Undir(:, 1:2); % Creates the matrix for the input

% Goes into the directory required
cd Computation_Module/SourceCode/Gopalan;

fid = fopen('Gopalan.txt', 'w'); % Creates text file
fprintf(fid, '%g\t%g\n', temp'); % Places data within file
fclose(fid);

% Runs command to find approximate number of communities
system(sprintf('./src/bin/svinet -file Gopalan.txt -n %g -k %g -label temp -findk'...
    , numnodes, numnodes));
% Adds the newly generated fold to the matlab path
addpath(genpath(pwd));

% Opens up the communities text file
fid = fopen('communities.txt', 'r');
test = textscan(fid, '%s', 'Delimiter', '\n'); % Imports the data in
fclose(fid);
system(sprintf('rm -rf n%g-k%g-temp-findk', numnodes, numnodes)); % Removes the folder

numComms = size(test{1}, 1); % Simply finds the number of rows for communities

% Runs the entire command to use the algorithm on the network
system(sprintf(...
    './src/bin/svinet -file Gopalan.txt -n %g -k %g -label Gopalan -link-sampling -max-iterations %g'...
    , numnodes, numComms, maxIter));

% Adds the newly generated folder to the matlab path
addpath(genpath(pwd));

system('rm Gopalan.txt'); % Deletes the network textfile (not needed any more)

% Creates the folder name of the output of Gopalan
foldername = sprintf('n%g-k%g-Gopalan-linksampling', numnodes, numComms);

% Moves the communities text file into the running directory
system(sprintf('mv %s/communities.txt ../../../Gopalan.txt', foldername));

% Deletes the old folder
system(sprintf('rm -rf %s', foldername));

% Goes into the directory required
cd ../../../;

fprintf('\n'); % Aesthetic value.
