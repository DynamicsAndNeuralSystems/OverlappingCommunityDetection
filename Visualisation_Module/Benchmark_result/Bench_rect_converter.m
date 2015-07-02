function[CommMat, I] = Bench_rect_converter(~, ~)
% Converts the format of the output of the benchmark test to be the same as
% other formats
% Brandon Lam, 10-07-2014
% Edits:
% 23-07-2014 - Changed the output into a matrix

%% Opens file, reads it.
fid = fopen('community.dat');
Commdata = textscan(fid, '%f\t%f %f %f %f %f %f %f %f %f %f');
fclose(fid); % This imports the dat file into MATLAB
BenchComm = cell(length(Commdata{1}), 1); % Makes a cell array

%% This for loop just creates a cell array in the same format as NodeLabels
for i = 1:size(Commdata{1}) % For the number of nodes
    n = 2;
    try while Commdata{n}(i,1) ~= 0 % While the next column isn't 0
              BenchComm{i} = [BenchComm{i}, Commdata{n}(i,1)]; % add it in
              n = n + 1;
        end
    catch % This is just so if it fails it keeps going because some rows don't have 0s
    end
    
    % Deletes all NaN from the cells
    BenchComm{i} = BenchComm{i}(~any(arrayfun(@isnan,BenchComm{i}), 3));
end

%% Changes the communities and orders them in terms of size
BenchComm = comm_replace(BenchComm);

%% Sorts nodes in community size
[~, I] = comm_sort(BenchComm);

%% Turns it into a matrix for rectangles to use
CommMat{1} = NodeLabel2Mat(BenchComm);