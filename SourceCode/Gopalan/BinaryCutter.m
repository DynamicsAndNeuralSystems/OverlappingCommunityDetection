% Cuts out links according to their weakness in the network
% Brandon Lam, 21-07-2014

% Takes input
prec = input('Please insert precision:\n');

load subject1.txt % Loads the input file

% Gets rid of all links weaker than the percentage of precision
net_mat = subject1(subject1(:,3)>quantile(subject1(:,3),prec), :);

% Makes a new file
file = ['network_' num2str(prec) '.txt'];
fid = fopen(file, 'wt');
fprintf(fid, '%d\t%d\n', net_mat(:,1:2)'); % Saves the output into this file
fclose(fid);