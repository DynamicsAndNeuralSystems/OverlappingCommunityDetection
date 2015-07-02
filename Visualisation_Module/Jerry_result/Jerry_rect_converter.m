function[CommMat] = Jerry_rect_converter(Bin, ~)

% Converts the output of the Jerry method into a form the rectangle
% visualiser can use
% Brandon Lam, 23-07-2014

for i = Bin
    % loads the file
    file = ['Jerry' num2str(i) '.mat'];
    load(file);
    
    % Sorts it to community size
    [Jerry] = comm_replace(Jerry);
    
    % Converts the nodelabels format to a matrix
    CommMat{Bin == i} = NodeLabel2Mat(Jerry);
end