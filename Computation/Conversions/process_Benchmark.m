function Output = process_Benchmark(benchFilename, numNodes)
% Function that processes the community data of the benchmark

% Reads the file
fid = fopen(benchFilename, 'r');
temp = textscan(fid, '%d', 'Delimiter', {' ', '\n'});
fclose(fid);

temp = temp{1}; % Easier access

zero_place = temp == 0;
% This code essentially cuts out every single 0, as well as the node
% numbers from the data, and then finds the maximum to find the number of
% communities
numComms = max(temp(logical(~zero_place.*~zero_place([end, 1:end-1]))));

% Finds where the node numbers are, used to therefore generate the final
% matrix for the communities
NodePlace = [find(zero_place([end, 1:end-1]) == 1); length(temp)+1];

% Preallocating for the calculations
Output = zeros(numNodes, numComms);

for i = 1:numNodes
    % Places values in the right place
    Output(i, temp(NodePlace(i)+1:NodePlace(i+1)-2)) =...
        1/length(NodePlace(i)+1:NodePlace(i+1)-2);
end

end
