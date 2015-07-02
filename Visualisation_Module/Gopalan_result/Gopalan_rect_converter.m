function[CommMat] = Gopalan_rect_converter(thresh, ~)

% Function that converts the output of the Gopalan algorithm to a format
% that the rectangle function will take and show in the visualiser.
% Brandon Lam, 23-07-2014

for i = thresh % Given a threshold range
    file = ['network_' num2str(i) '_groups.txt']; % Gets the file path
    groups = load(file); % Loads the entire matrix
    [~, I] = sort(groups(:, 2)); % Sorts the second column
    temp = groups(I, 3:end); % Fixes everything up to be a good matrix
    % Also gets rid of all unneeded numbers
    
    [~, Ind] = sort(sum(temp)); % Finds the order of communities, and changes the order around
    
    CommMat{thresh == i} = temp(:, Ind); % Adds it to the big cell
end