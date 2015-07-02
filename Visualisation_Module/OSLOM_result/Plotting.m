% Brandon Lam, 07-07-2014
% Takes data from OSLOM analysed data and converts it to three plots
% Plot 1 - Tolerance against number of communities
% Plot 2 - Tolerance against number of overlapping nodes
% Plot 3 - Tolerance against number of homeless nodes

clear all; clc; close all;

% These will have to be changed if this is turned into a function
NumSubject = 1; % Put the number of subjects
Tolerances = 0.1:0.1:1; % Range of tolerances

% Preallocating for speed
NumCommMat = zeros(NumSubject, length(Tolerances));
NumOverlap = zeros(NumSubject, length(Tolerances));
NumHomeles = zeros(NumSubject, length(Tolerances));

% This loop goes through all subjects, taking their data and making
% matrices out of the data
for i = 1:NumSubject
    % Dir = ['subject' num2str(i)]; % This is the directory where it takes the txt files from
    for tol = Tolerances
        File = ['subject' num2str(i) '_data_tol_' num2str(tol) '.txt']; % This is the file name
        % Path = [Dir '/' File]; % The path is therefore constructed
        [NodeLabels, NumComms] = Converter(File); % Calls upon my conversion function
        
        % Creation of the Number of communities matrix
        NumCommMat(1, Tolerances == tol) = NumComms;
        
        % Creation of the overlapping node communities
        NumOverlap(1, Tolerances == tol) = sum(cellfun('length', NodeLabels) > 1);
        
        % Creation of the homeless nodes matrix
        NumHomeles(1, Tolerances == tol) = sum(cellfun('length', NodeLabels) < 1);
    end
end

figure('units','normalized','outerposition',[0 0 1 1]) % Makes a full screen figure
subplot(2, 2, 1); % Creates 4 subplots to use for the graphs
% Number of communities plot
plot(Tolerances, NumCommMat, 'k');
title('Average Number of Communities');
xlabel('Tolerances'); ylabel('Number of Communities');

% Number of overlapping nodes plot
subplot(222);
plot(Tolerances, NumOverlap, 'r')
title('Average Number of Overlapping Nodes');
xlabel('Tolerances'); ylabel('Number of Overlapping Nodes');

% Number of homeless nodes plot
subplot(223);
plot(Tolerances, NumHomeles, 'r')
title('Average Number of Homeless Nodes');
xlabel('Tolerances'); ylabel('Number of Homeless Nodes');