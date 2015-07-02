function[CommMat] = OSLOM_rect_converter(Tolerances, ~)
% Runs reiteratively to convert the OSLOM outputs into Nodelabels which is
% useful for comparing algorithms
% Brandon Lam, 10-07-2014
% Edits:
% 23-07-2014 - Changes the output into matrices within a cell

% Tolerances = 0.1:0.1:1; % Range of tolerances

% This loop runs the conversion function for each tolerance level
for i = Tolerances
    file = ['subject1_data_tol_' num2str(i) '.txt']; % Defines the name
    [NodeLabels, ~] = Converter(file); % Calls function to obtain from data
    
    NodeLabels = comm_replace(NodeLabels); % Calls function to change community values
    
    temp = NodeLabel2Mat(NodeLabels); % temporarily creates a matrix to put in the big one
            
    CommMat{Tolerances == i} = temp; % Adds it to the matrix
%     NumComms(1, Tolerances == i) = NumCommsTemp; % Shows the number of communities
end