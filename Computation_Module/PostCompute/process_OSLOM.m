function [Output] = process_OSLOM(Tol, numnodes)
% Function that processes the output of the OSLOM method

Output = cell(1, length(Tol)); % Creates the cells for the output

for t = Tol
    % Defines the filename for program to read
    filename = sprintf('OSLOM_tol_%g.txt', t);
    
    fid = fopen(filename); % Opens file
    % Scans the text, creating a list, with 0s where it is the next module
    A = textscan(fid, '%d', 'CommentStyle', '#', 'Delimiter', {' ', '\n'});
    A = A{1}; % Sets it as the matrix, gets rid of the cell
    
    % Closes and deletes the textfile
    fclose(fid); 
    command = ['rm ' filename]; % Creates the command to delete the file
    system(command); % Deletes
    
    NumComms = sum(A==0); % Finds the number of communities
    
    CommMat = zeros(numnodes, NumComms); % Creates the matrix for each tolerance
    
    counter = 1; % Counter for communities
    
    for i = 1:length(A) % For the entire A list
        if ~A(i) == 0 % If it isn't 0
            CommMat(A(i), counter) = 1; % Places a 1 in the community matrix
        else
            counter = counter + 1; % Else increases the counter
        end
    end
    
    for i = 1:numnodes
        % Equally distributes the communities
        CommMat(i, :) = CommMat(i, :)/sum(CommMat(i, :));
    end
    
    Output{Tol == t} = CommMat; % Sets it as the output
end