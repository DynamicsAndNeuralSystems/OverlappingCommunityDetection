function[CommMat] = Link_rect_converter(thresh, numnodes)
% Converts the output of the link communities code into a form that my
% visualiser can use
% Brandon Lam, 16-07-2014
% Edits:
% 23-07-2014 - Changed the output into matrices within a cell

for i = thresh
    file = ['Comm2nodes_' num2str(i) '.txt'];
    fid = fopen(file); % Opens the text file
    S = textscan(fid, '%s', 'Delimiter','\n'); % Reads the string lines
    fclose(fid); % Closes textfile
    
    NumLinkComms = length(S{1}); % Finds the number of communities
    % % Tells the user information
    % display(['We found ' num2str(NumLineComms) ' communities'])
    
    % Preallocates cells
    D = cell(NumLinkComms,1);
    
    % Loop that deletes the first value of all lines, and puts it into cells
    for j = 1:NumLinkComms
        SS = textscan(S{1}{j},'%f'); % Reads each line
        D{j} = SS{1}(2:end); % New cell is missing the first value
        %     % Tells the user how many nodes are in each community
        %     fprintf(1,'There are %u nodes in community %u\n',length(D{i}),i);
    end
    
    % % Sorts the communities in terms of sizes
    % [~,ix] = sort(cellfun(@length,D),'descend');
    % D = D(ix); % Changes D into the new D
    
    % Defines the new NodeLinkLabels
    temp = arrayfun(@(x)find(cellfun(@(y)ismember(x,y),D))',1:numnodes,'UniformOutput',0)';
    temp = comm_replace(temp); % Calls function to change community values
    
    % Calls function to turn it into a matrix
    Mat = NodeLabel2Mat(temp);
    
    CommMat{thresh==i} = Mat; % Adds that matrix to the whole cell
end