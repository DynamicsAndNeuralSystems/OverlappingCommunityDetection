function [Output] = process_Link(numnodes)
% Processes the textfile output of the Gopalan Method.

% Scans the numbers within it
temp = importdata('Link.txt', '\n');

% Deletes the file
system('rm Link.txt');

numComms = length(temp); % Finds the number of communities

% Preallocates cells
Comms = cell(numComms,1);

% Loop that deletes the first value of all lines, and puts it into cells
for j = 1:numComms
    temp2 = textscan(temp{j},'%f'); % Reads each line
    Comms{j} = temp2{1}(2:end); % New cell is missing the first value
end

% Defines the new NodeLinkLabels
temp3 = arrayfun(@(x)find(cellfun(@(y)ismember(x,y),Comms))',1:numnodes,'UniformOutput',0)';

% Preallocates the matrix for the output
Output = zeros(numnodes, numComms);

for i = 1:numnodes
    val = 1/length(temp3{i}); % Calculates the value to put in (equal amounts)
    Output(i, temp3{i}) = val; % Finally places those values 
end