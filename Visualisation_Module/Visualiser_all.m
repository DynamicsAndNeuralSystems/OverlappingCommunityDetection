%% Visualiser to compare the different methods of overlapping network
% detection. A benchmark is used as the reference.
% Brandon Lam, 10-07-2014

clear all; close all; clc; run ./startup.m

%% Variables that can be changed
Width = 4; % Only integers please
Methods = {'Bench', 'OSLOM', 'Link'};
figure('color', 'w'); cmp = colormap(); % Put your favourite colour map here

% Various inputs
Gopalan = [0 0.4]; % Thresholds done on network
OSLOM = 0.1:0.1:1; % Thresholds within algorithm
Link = [0 0.1 0.4 0.6 0.8 0.9]; % Thresholds done on network
Jerry = [0 1]; % Unweighted or Weighted
Shen = 2:10; % Clique sizes of K

%% Obtains the adjacency matrix
[numnodes, Adj_Mat] = Ntwk2AdjMat; % Using the benchmark as the network

%% Read community.dat to get the communities made from the benchmarking process
[Sort_Bench_Mat, I] = Bench_rect_converter;

%% Sorting in terms of communities
Sort_Adj_Mat = Adj_Mat(I, I); % Converts the original Adjacency matrix to the community ordered one
Sort_Adj_Mat = Sort_Adj_Mat/max(max(Sort_Adj_Mat)); % Normalise adjacency matrix
Sort_Adj_Mat(Sort_Adj_Mat == 0) = NaN; % Sets the 0s to NaNs so the program will set them black

StartPos = size(Sort_Adj_Mat, 2) + 0.5; % Sets the starting x position for the rectangles

%% Plotting
for Case = Methods % Loops for each method
    switch Case{1} % Defines variables for the functions that require them
        case 'Gopalan';     Var = Gopalan; % Thresholds of network
        case 'OSLOM';       Var = OSLOM; % Thresholds within function
        case 'Link';        Var = Link; % Thresholds of network
        case 'Jerry';       Var = Jerry; % Unweighted and Weighted
        case 'Shen';        Var = Shen; % Different k values
        otherwise;          Var = 0; % Random thing, can be anything
    end
    % Adds to the size of the matrix
    Sort_Adj_Mat = [Sort_Adj_Mat NaN(numnodes, size(Var, 2)*Width)];
end

Plot = imagesc(Sort_Adj_Mat); hold on; % Plots the image with colour
set(gca, 'color', [0 0 0], 'CLim', [0 1]); % Sets background colour to black
set(Plot, 'alphadata',~isnan(Sort_Adj_Mat)); % Turns all NaN values into transparent colours
colorbar('location', 'eastoutside'); % Shows a colour bar
title('Comparisons of algorithms');
xlabel('Algorithms');

%% Running converters for all methods specified
LabelsStr = cell(0, 1);
for Case = Methods % Loops for each method
    Conv = [Case{1} '_rect_converter']; % Defines the name of the function
    Func = str2func(Conv); % Creates a function handle for later use
    switch Case{1} % Defines variables for the functions that require them
        case 'Gopalan';     Var = Gopalan; % Thresholds of network
        case 'OSLOM';       Var = OSLOM; % Thresholds within function
        case 'Link';        Var = Link; % Thresholds of network
        case 'Jerry';       Var = Jerry; % Unweighted and Weighted
        case 'Shen';        Var = Shen; % Different k values
        otherwise;          Var = 0; % Random thing, can be anything
    end
    
    CommMat = Func(Var, numnodes); % Calls the specific function to use
    for num = 1:length(CommMat) % Loops for all outputs (since there can be multiple matrices)
        % Plots a black line for every different method/setting
        plot([StartPos, StartPos], [0, numnodes + 1], 'k');
        CommMat{num} = CommMat{num}(I, :); % Sorts the matrix first
        for i = 1:numnodes % Loops for all nodes
            TempPos = StartPos; % Resets the temporary starting position
            for j = 1:size(CommMat{num}, 2) % Loops for all communities
                if CommMat{num}(i, j) > 0
                    % Creates a rectangle for each community they exist in,
                    % with weightings as well. The colour ischosen from a
                    % 64 size colour map, by multiplying 64 by the
                    % fraction, then rounding to the next integer.
                    rectangle('Position', [TempPos, i - 0.5, Width*CommMat{num}(i, j), 1], ...
                        'FaceColor', cmp(ceil(64*(j/size(CommMat{num}, 2))), :), ...
                        'EdgeColor', 'None');
                    % Changes the temporary position to allow for next rectangle
                    TempPos = TempPos + Width*CommMat{num}(i, j);
                end
            end
        end
        LabelsStr{end+1} = sprintf('%s%g',Case{1}, Var(num)); % Adds this to the labels
        StartPos = StartPos + Width; % Changes the starting position, for next algorithm/output
        if size(CommMat{num}, 2) > 64
            % Warns the user of possible problems with colours
            warning('%s %g has too many communities to be accurately displayed!', Case{1}, Var(num));
        end
    end
end

%% Labelling
XTickLabel = size(Sort_Adj_Mat, 1)+0.5; % Kickstarts the for loop
for i = 1: size(LabelsStr, 2)
    XTickLabel = [XTickLabel, XTickLabel(end)+Width]; % Finds all the labels needed
end
set(gca, 'XTick', XTickLabel+Width/2); % Sets the axis ticks to these values
set(gca, 'XTickLabel', LabelsStr, 'FontSize', 7); % Shows the labels, makes them smaller
