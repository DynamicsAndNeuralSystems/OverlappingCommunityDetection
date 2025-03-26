function Output = process_OSLOM(Tol,numNodes, oslomResPath)
% Function that processes the output of the OSLOM method

Output = cell(1, length(Tol)); % Creates the cells for the output

for t = Tol
    % Defines the fileName for program to read
    fileName = sprintf('%s/OSLOM_tol_%g.txt', oslomResPath, t);

    fid = fopen(fileName); % Opens file
    % Scans the text, creating a list, with 0s where it is the next module
    A = textscan(fid, '%d', 'CommentStyle', '#', 'Delimiter', {' ', '\n'});
    A = A{1}; % Sets it as the matrix, gets rid of the cell

    % Closes and deletes the textfile
    fclose(fid);
    command = sprintf('rm %s',fileName);
    system(command); % Deletes the file

    numComms = sum(A==0);

    % Creates the matrix for each tolerance:
    commMat = zeros(numNodes, numComms);

    % Counter for communities
    counter = 1;

    for i = 1:length(A) % For the entire A list
        if ~A(i) == 0 % If it isn't 0
            commMat(A(i), counter) = 1; % Places a 1 in the community matrix
        else
            counter = counter + 1; % Else increases the counter
        end
    end

    for i = 1:numNodes
        % Equally distributes the communities
        commMat(i, :) = commMat(i, :)/sum(commMat(i, :));
    end

    Output{Tol == t} = commMat; % Sets it as the output
end
