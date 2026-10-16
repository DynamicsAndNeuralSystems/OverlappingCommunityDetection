function [fig_h, ax_h] = OLCD_Visualization(Methods, resultsFile, cmap)
% Visualises results of the Computation Module: the (sorted) network adjacency
% matrix, followed by one column per method output showing each node's
% community memberships as colored bars.
%-------------------------------------------------------------------------------
% Inputs:
%   Methods:     cell array (or string array/char) of method names to plot.
%                Each entry matches either the method's .Name (e.g., 'NNMF'
%                plots all NNMF_* outputs) or an exact field name of Final
%                (e.g., 'NNMF_20'). Names matching nothing give a warning.
%                Default (omitted/empty): all method outputs in the file
%                (every field except Date and Network, including Benchmark).
%   resultsFile: results file containing the structure Final
%                (default: 'Computation_Result.mat', must be on the path)
%   cmap:        colormap: a name ('jet'), function handle (@parula), or an
%                N x 3 matrix (default: 'parula')
%
% Outputs:
%   fig_h: handle to the figure
%   ax_h:  handle to the axes
%-------------------------------------------------------------------------------

if nargin < 2 || isempty(resultsFile)
    resultsFile = 'Computation_Result.mat';
end
if nargin < 3 || isempty(cmap)
    cmap = 'parula';
end

%% Load data
loaded = load(resultsFile,'Final');
Final = loaded.Final;
fprintf('Loading data created on %s!\n', char(string(Final.Date)));
numNodes = size(Final.Network, 1);

%% Resolve colormap
if ischar(cmap) || isstring(cmap)
    cmp = feval(char(cmap), 256);
elseif isa(cmap, 'function_handle')
    cmp = cmap(256);
else
    cmp = cmap;
end
if size(cmp, 2) ~= 3
    error('cmap must be a colormap name, function handle, or N x 3 matrix');
end
numColors = size(cmp, 1);

%% Finding the methods within Final
allFields = fieldnames(Final);
allFields = allFields(~ismember(allFields, {'Date','Network'}));
% Keep only method outputs (structs with Result):
isMethod = cellfun(@(f) isstruct(Final.(f)) && isfield(Final.(f),'Result'), allFields);
allFields = allFields(isMethod);

if nargin < 1 || isempty(Methods)
    methodNames = allFields';
else
    Methods = cellstr(Methods); % Allows char/string input
    methodNames = cell(0); % Method names used
    for Vis_name = Methods(:)'
        matched = false;
        for i = 1:length(allFields)
            if strcmp(Final.(allFields{i}).Name, Vis_name{1}) || strcmp(allFields{i}, Vis_name{1})
                matched = true;
                if ~ismember(allFields{i}, methodNames)
                    methodNames{end+1} = allFields{i}; %#ok<AGROW>
                end
            end
        end
        if ~matched
            warning('No results found matching ''%s''', Vis_name{1});
        end
    end
end
numMethods = length(methodNames);
disp('Found all results pertaining to the input.');

%% Set parameters
Width = 6; % Integers, defines width of each column
fig_h = figure('color', 'w');
colormap(fig_h, cmp);
% Sensible figure size, kept within the screen:
scr = get(groot, 'ScreenSize');
figW = min(1536, scr(3)); figH = min(703, scr(4)-50);
fig_h.Position = [scr(1) + (scr(3)-figW)/2, scr(2) + (scr(4)-figH)/2, figW, figH];

%% Benchmark sorting of nodes
if isfield(Final, 'Benchmark') && isfield(Final.Benchmark, 'Result')
    sortBy = Final.Benchmark.Result;
    disp('Nodes have been sorted into communities (benchmark)');
elseif numMethods > 0
    sortBy = Final.(methodNames{1}).Result;
    fprintf('No benchmark found: nodes sorted by %s\n', methodNames{1});
else
    sortBy = zeros(numNodes, 0);
    disp('No benchmark found: nodes in natural order');
end
I = Node_Sorter(sortBy); % Sorts the nodes into communities

%-------------------------------------------------------------------------------
%% Initial matrix view
%-------------------------------------------------------------------------------
full_Matrix = Final.Network/max(max(Final.Network)); % Normalises the data for plotting
full_Matrix = full_Matrix(I, I); % Moves the nodes to their sorted locations
full_Matrix(full_Matrix == 0) = NaN; % Turns the 0s to NaNs for easy plotting
disp('Plotted network matrix.');
StartPos = size(full_Matrix, 2) + 0.5; % Starting position of x position

% Add blank space for the plot, for each method
full_Matrix = [full_Matrix, NaN(size(full_Matrix, 1), Width*numMethods)];

%-------------------------------------------------------------------------------
%% Figure options
%-------------------------------------------------------------------------------
Plot = imagesc(full_Matrix); hold on
ax_h = gca;
set(ax_h, 'color', [0 0 0], 'CLim', [0 1]); % Sets background colour to black
set(Plot, 'alphadata',~isnan(full_Matrix)); % Turns all NaN values into transparent colours
colorbar('location', 'eastoutside'); % Shows a colour bar

%-------------------------------------------------------------------------------
%% Plotting rectangles and lines
%-------------------------------------------------------------------------------
for name = methodNames
    % Plots a starting line, so that algorithms can be separated
    plot([StartPos, StartPos], [0, numNodes + 1], 'k');

    % Temporary community matrix for each algorithm
    CommMat = Final.(name{1}).Result(I, :);

    % Reorders nodes to fit colours to largest to smallest (NaN treated as 0)
    CommMat = Node_Reorder(CommMat);
    numComms = size(CommMat, 2);

    % One rectangle per nonzero entry, drawn as a single patch object.
    % Each rectangle shows the strength of connection to a community.
    rightEdge = StartPos + Width*cumsum(CommMat, 2);
    leftEdge = rightEdge - Width*CommMat;
    [yInd, xInd] = find(CommMat ~= 0);
    lin = sub2ind(size(CommMat), yInd, xInd);
    nRect = length(lin);
    if nRect > 0
        xL = leftEdge(lin); xR = rightEdge(lin);
        Vx = [xL, xR, xR, xL]';
        Vy = [yInd-0.5, yInd-0.5, yInd+0.5, yInd+0.5]';
        Vertices = [Vx(:), Vy(:)];
        Faces = reshape(1:4*nRect, 4, nRect)';
        colorInd = min(max(ceil(numColors*(xInd/numComms)), 1), numColors);
        patch('Faces', Faces, 'Vertices', Vertices, ...
            'FaceVertexCData', cmp(colorInd, :), 'FaceColor', 'flat', ...
            'EdgeColor', 'none');
    end

    if numComms > numColors % Alerts the user if there are more communities than colors
        warning('%s method has more communities than colors in the colormap, so colors are repeated!', name{1});
    end
    fprintf('%s method has been plotted.\n', name{1});
    StartPos = StartPos + Width; % Changes the starting position, for next algorithm/output
end

%% Labelling
XTick = numNodes + 0.5 + Width*((1:numMethods) - 0.5); % Center of each method's column
set(ax_h, 'XTick', XTick); % Sets the axis ticks to these values
labels = strrep(methodNames, '_', '\_');
set(ax_h, 'XTickLabel', labels, 'FontSize', 10); % Shows the labels, makes them smaller
set(ax_h, 'XTickLabelRotation', 45); % Rotates the labels to 45 degrees
set(fig_h, 'InvertHardCopy', 'off'); % Fixes white background issue

end
