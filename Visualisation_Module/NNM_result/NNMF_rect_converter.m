function[CommMat] = NNMF_rect_converter(~, ~)

% Converts the output of the NNMF algorithm into something the visualiser
% can use
% Brandon Lam, 23-07-2014

% Loads the .matfile
load('NNMF.mat');

[~, Ind] = sort(sum(NNMF)); % Finds the order of communities, and changes the order around

CommMat{1} = NNMF(:, Ind); % Adds it to the big cell
