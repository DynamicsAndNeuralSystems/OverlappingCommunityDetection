function [Output] = process_infomap(numnodes)
% Function that processes the output of the Infomap method into the matrix

cd output
Output = cell(1, 1);

% Defines the filename for program to read
filename = 'networki.clu';
fid = fopen(filename); % Opens file
% Scans the text, creating a list
A = textscan(fid, '%d %d %d', 'CommentStyle', '#', 'Delimiter', {' ', ' ',' ', '\n'});
B = A{1}; % Sets it as the matrix, gets rid of the cell
C = A{2};
A = [B,C];

% Closes and deletes the textfile
fclose(fid); 
%     command = ['rm ' filename]; % Creates the command to delete the file
%     system(command); % Deletes

NumComms = max(A(:,2)); % Finds the number of communities
CommMat = zeros(numnodes, NumComms); % Creates the matrix for each tolerance

for j = 1:length(A) % For the entire A list
    CommMat(A(j,1),A(j,2))=1; 
end

for j = 1:numnodes
    % Equally distributes the communities
    CommMat(j, :) = CommMat(j, :)/sum(CommMat(j, :));
end

Output = CommMat; % Sets it as the output

cd ../
system('rm -r output')
cd ../../
end