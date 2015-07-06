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
maxIter = 1000; % Maximum iterations done within the algorithm
prec_Gopalan = [0 0.1 0.2 0.3]; % Percentage of lowest links to cut before computation

% ITERATIONS OF PRECISION
Gopalan_final = cell(0); % Preallocates the cell for output

for c = prec_Gopalan
    % Gets rid of all links weaker than the percentage of precision
    Input = Undir(Undir(:,3)>quantile(Undir(:,3),c), :);
    
    % Runs the Gopalan method
    run_Gopalan(Input, numnodes, maxIter);
    
    % Processes the data
    [Gopalan_output] = process_Gopalan(numnodes);
    
    % Places the data within the final output
    Gopalan_final{1, end + 1} = Gopalan_output;
end

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

% Processes the output
[Jerry_final] = process_Jerry(Jerry, numnodes); 

%% Link
% Input is undirected list
% ------------------------------PARAMETERS---------------------------------
prec_Link = [0 0.01 0.1]; % Percentage of lowest links to cut before computation

% ITERATIONS OF PRECISION
Link_final = cell(0); % Preallocates the cell for output

for c = prec_Link
    % Gets rid of all links weaker than the percentage of precision
    Input = Undir(Undir(:,3)>quantile(Undir(:,3),c), :);
    
    % Runs the link communities method
    run_Link(Input);
    
    % Processes the textfile outputs
    Link_output = process_Link(numnodes);
    
    % Places the data within the final output
    Link_final{1, end + 1} = Link_output;
end

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

%% Constructing the structure
Computation = struct([]);












