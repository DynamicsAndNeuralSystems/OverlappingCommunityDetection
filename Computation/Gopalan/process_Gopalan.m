function [Output] = process_Gopalan(numnodes)
% Processes the textfile output of the Gopalan Method.

% Opens the text file
fid = fopen('Gopalan.txt', 'r');

% Scans the numbers within it
temp = textscan(fid, '%d', 'Delimiter', {' ', '\n'});

% Closes and deletes the file
fclose(fid);
system('rm Gopalan.txt');

temp = temp{1}; % Easier to access

counter = 1; % Counter for communities
numComms = sum(temp == 0); % Counts the number of communities

Output = zeros(numnodes, numComms); % Preallocates the matrix

% For every number in temp
for i = 1:length(temp)
    if temp(i) == 0
        counter = counter + 1; % Adds one to the counter if it encounters a 0
    else
        Output(temp(i), counter) = 1; % Sets a 1 in a specific place
    end
end

for x = 1:numnodes
    % divides the values according to how many 1s are in it
    Output(x, :) = Output(x, :)/sum(Output(x, :)); 
end