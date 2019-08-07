function [Plot] = Visualisebrain(network)
%% Visualiser code
% Brandon Lam, 07-07-2015
% Visualisation.m edited to visualise a network from its adjmat

% clc; clear all; close all;

%% Loading data
Final.Network = network;
numnodes = size(Final.Network, 1);

%% Parameters
fig_h = figure('color', 'w');
cmp = colormap(jet); % Put your favourite colour map here
% fig_h.Position = [100, 100, 1700, 1000];
%fig_h.Position = [1, 26, 1536, 703];

%% Initial matrix view
full_Matrix = Final.Network/max(max(Final.Network));
full_Matrix(full_Matrix == 0) = NaN; % Turns the 0s to NaNs for easy plotting
disp('Plotted network matrix.');

%% Figure options
subplot(2,2,[1 2]);
Plot = imagesc(full_Matrix); hold on
axis image;
set(gca, 'color', [0 0 0], 'CLim', [0 1]); % Sets background colour to black
set(Plot, 'alphadata',~isnan(full_Matrix)); % Turns all NaN values into transparent colours
colorbar('location', 'eastoutside'); % Shows a colour bar
%% Labelling
title('A','FontSize', 15);
XTickLabel = size(full_Matrix, 1)+0.5; % Kickstarts the for loop
set(gcf, 'InvertHardCopy', 'off'); % Fixes white background issue

end