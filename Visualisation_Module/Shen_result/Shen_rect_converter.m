function[CommMat] = Shen_rect_converter(Range, ~)
% Takes the data from the Shen data and converts it to a node labels format
% Brandon Lam, 18-07-2014
% Edits:
% 23-07-2014 - changes the output into a matrix

% Range = 2:10; % Put the range of your k values
Thr = 1.3; % Put the threshold you want to see

for i = Range
    % Defines the file path
    file = ['Weight_Thresh_' num2str(Thr) '/Shen_BenchMark_k_' num2str(i) '.mat'];
    
    % Loads the file into MATLAB
    thisLoad = load(file);
    
    % Creates a nodelabels, and replaces the communities with sorted size
    NodeLabels = comm_replace(thisLoad.NodeLabels);
    
    % Turns this into a matrix
    temp = NodeLabel2Mat(NodeLabels);
    
    CommMat{Range==i} = temp;
end
