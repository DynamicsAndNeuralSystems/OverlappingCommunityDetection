function [DirList] = Mat2Direct(Mat, numnodes)
% Converts a matrix to a directed list

DirList = zeros(sum(sum(Mat ~= 0)), 3);
counter = 1; % Counter to put values into the directed list

for x = 1:numnodes
    for y = 1:numnodes
        if Mat(x, y) ~= 0
            % Places the values into the directed list
            DirList(counter, :) = [x, y, Mat(x, y)];
            counter = counter + 1; % Adds 1 to the counter for the next values
        end
    end
end
