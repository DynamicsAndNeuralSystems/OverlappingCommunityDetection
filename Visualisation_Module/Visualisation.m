%% Visualiser code
% Program that takes input from the computational module and outputs a
% visualisation of the results
% Brandon Lam, 07-07-2015

clc; clear all; close all;

%% Parameters
Width = 4; % Integers, defines width of each column
fig_h = figure('color', 'w'); cmp = colormap(jet); % Put your favourite colour map here
fig_h.Position = [174, 133, 1686, 986];
Methods = {'Benchmark', 'OSLOM', 'Shen', 'Jerry'}; % Put your favourite methods here!
MethodName = cell(0); % Method names used


%% Loading data
load('Computation_Result');
fprintf('Loading data created on %s!\n', Final.Date);
numnodes = size(Final.Network, 1);

%% Benchmark sorting of nodes
[I] = Node_Sorter(Final.Benchmark.Result); % Sorts the nodes into communities
disp('Nodes have been sorted into communities');

%% Initial matrix view
full_Matrix = Final.Network/max(max(Final.Network)); % Normalises the data for plotting
full_Matrix = full_Matrix(I, I); % Moves the nodes to their sorted locations
full_Matrix(full_Matrix == 0) = NaN; % Turns the 0s to NaNs for easy plotting
disp('Plotted network matrix.');
StartPos = size(full_Matrix, 2) + 0.5; % Starting position of x position


%% Finding the methods within Final
for Vis_name = Methods
    % Temporary data dump of names
    temp = fieldnames(Final);
    temp = temp(3:end); % Gets rid of 'Date' and 'Network' - not needed
    
    for i = 1:length(temp)
        
        if strcmp(Final.(temp{i}).Name, Vis_name{1}) % - If the names match
            % Add blank space for the plot
            full_Matrix = [full_Matrix NaN(size(full_Matrix, 1), Width)];
            % Saves the full name of the method within the "Final"
            % structure
            MethodName{end+1} = temp{i};
        end
    end
end

disp('Found all results pertaining to the input.');

%% Figure options
Plot = imagesc(full_Matrix); hold on
set(gca, 'color', [0 0 0], 'CLim', [0 1]); % Sets background colour to black
set(Plot, 'alphadata',~isnan(full_Matrix)); % Turns all NaN values into transparent colours
colorbar('location', 'eastoutside'); % Shows a colour bar


%% Plotting rectangles and lines
for name = MethodName
    % Plots a starting line, so that algorithms can be separated
    plot([StartPos, StartPos], [0, numnodes + 1], 'k');
    CommMat = Final.(name{1}).Result(I, :);
    
    for y = 1:numnodes
        
        tempPos = StartPos; % Resets the temporary position
        
        for x = 1:size(CommMat, 2)
            if isnan(CommMat(y, x))
                break
            end
            % Creates rectangles to show the strength of connection to a
            % community (This cycles through every value to create the
            % rectangles
            if ~CommMat(y, x) == 0 % Skips the 0s. Speeds up the process quite a bit
                rectangle('Position', [tempPos, y - 0.5, Width*CommMat(y, x), 1], ...
                    'FaceColor', cmp(ceil(64*(x/size(CommMat, 2))), :), ...
                    'EdgeColor', 'None');
            end
            
            tempPos = tempPos + Width*CommMat(y, x); % Sets the next position
        end
    end
    
    if size(CommMat, 2) > 64 % Alerts the user if there are too many communities
        warning('%s method has too many communities to be plotted accurately!', name{1});
    end
    fprintf('%s method has been plotted.\n', name{1});
    StartPos = StartPos + Width; % Changes the starting position, for next algorithm/output
end

%% Labelling
XTickLabel = size(full_Matrix, 1)+0.5; % Kickstarts the for loop
for i = 1: size(MethodName, 2)
    XTickLabel = [XTickLabel, XTickLabel(end)+Width]; % Finds all the labels needed
end
set(gca, 'XTick', XTickLabel+Width/2); % Sets the axis ticks to these values
set(gca, 'XTickLabel', MethodName, 'FontSize', 7); % Shows the labels, makes them smaller