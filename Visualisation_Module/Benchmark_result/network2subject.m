% Converts the network to a undirected network matrix
% Brandon Lam, 09-07-2014

load network.dat % Loads the file
% Gets rid of everything that is higher than the second node
A = network(network(:, 1)<network(:, 2), :);
% Opens a new file to write
fid = fopen('subject1.txt', 'wt');
fprintf(fid, '%d\t%d\t%.4f\n', A'); % Writes it in
fclose(fid);