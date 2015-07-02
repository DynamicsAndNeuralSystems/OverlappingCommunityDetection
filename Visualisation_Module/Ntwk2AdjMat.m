function [numnodes, Adj_Mat] = Ntwk2AdjMat
% m file that converts the network.dat file into an adjacency matrix
% Brandon Lam, 10-07-2014

load ./Benchmark_result/network.dat % Loads the network file into MATLAB
numnodes = max(network(:, 2)); % Calculates the number of nodes
Adj_Mat=zeros(numnodes,numnodes); % Preallocates the adjacency matrix

for i = 1:size(network) % For every entry in the network file
    % It adds the value into the matrix
    Adj_Mat(network(i, 1), network(i, 2)) = network(i, 3); 
end