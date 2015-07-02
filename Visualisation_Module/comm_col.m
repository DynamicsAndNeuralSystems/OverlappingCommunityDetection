function[comms] = comm_col(sorted_comm, numcomm)
% Takes in a community cell and the number of communities, and outputs a
% matrix that can be used to compare against others by adding colours.
% Brandon Lam, 11-07-2014

% Sets each label to double first (if multiple)
comms = cellfun(@(x)double(x(1)),sorted_comm);
% Over-writes labels for nodes with multiple membership:
comms(cellfun(@length,sorted_comm)>1) = max(numcomm)+1;
% Normalize to [0,1] interval
comms = (comms-min(comms))./(max(comms)-min(comms));
end