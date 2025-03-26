function [] = run_infomap(dirlist,sourceCodePath)
% Function that runs the link communities command

% Goes into the directory required
cd(fullfile(sourceCodePath, 'Infomap'));

fid = fopen('networki.txt', 'w'); % Creates text file
fprintf(fid, '%g\t%g\t%f\n', dirlist'); % Places data within file
fclose(fid);

system('rm -r output')
mkdir output

system('./Infomap networki.txt output --overlapping --clu ');
 system('rm networki.txt')
