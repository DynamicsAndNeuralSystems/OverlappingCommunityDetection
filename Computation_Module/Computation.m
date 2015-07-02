% Program that computes algorithms for a specific network
% Brandon Lam, 29-06-2015

clc; close all; clear all; startup;

%% Input the directed data list here
DirList = load('network.dat'); % Loads the directed data
numnodes = max(DirList(:, 1)); % Calculates the number of nodes in the system

%% Creates the input files
Mat = Direct2Matrix(DirList, numnodes); % Calls function to make matrix input
Undir = Direct2Undir(DirList); % Calls function to make undirected input

%% Clauset - uses weighted directed list as an input
% This method is relatively easy - it uses a directed list of connections
% as an input, but the output has no overlaps.

Clauset = run_Clauset(DirList); % Runs the Clauset method
Clauset_final = process_Clauset(Clauset, numnodes); % Processes the output

%% Gopalan
% Input is either undirected or directed list, but it will treat them as
% undirected
% ------------------------------PARAMETERS---------------------------------



%% Jerry
% Input is Adjacency matrix
% ------------------------------PARAMETERS---------------------------------
numIters = 120;
Threshold = 0.09;
IsBlind = 1;
ListeningRule = 'majority'; 
% ListeningRule = 'probabilistic'; 

% Runs the Jerry method
[Jerry] = run_Jerry(Mat, numIters, Threshold, IsBlind, ListeningRule, numnodes);
[Jerry_final] = process_Jerry(Jerry, numnodes); % Processes the output

%% Link



%% NNMF
% Input is Adjacency matrix
% ------------------------------PARAMETERS---------------------------------
thresh = [0 0.01 0.1]; % Thresholds, can be multiple thresholds.

NNMF = run_NNMF(Mat, thresh); % Runs the NNMF method

% This has no need for processing, as the output is the required matrix

%% OSLOM
% Input is a sparse matrix (undirected)
% ------------------------------PARAMETERS---------------------------------
numIter = 100; % Number of iterations within the algorithm
Tol = 0.1:0.1:1; % Range of tolerances
filepath = './Modules/Computation_Module/Computation/OSLOM'; % File path from module to OSLOM code
fileback = '../../../../'; % File path back

% Runs the algorithm for OSLOM
run_OSLOM(Undir, numIter, Tol, filepath, fileback);

% Processes the textfile outputs
[OSLOM_final] = process_OSLOM(Tol, numnodes);

%% Shen
% Input is Adjacency matrix
% ------------------------------PARAMETERS---------------------------------
cSize = 7; % Clique size
cwThresh = 1.3; % Clique weight threshold
lThresh = 0; % Threshold on links, 0 because we want to preserve data

% Runs the function for the Shen method
[Shen] = run_Shen(Mat, cSize, cwThresh, lThresh);

% Processes the data into the matrix
[Shen_final] = process_Shen(Shen, numnodes);

