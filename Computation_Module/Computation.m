% clc; close all; clear all; startup;
function Final = Computation(input,isBenchmark)
% Program that computes algorithms for a specific network
% Brandon Lam, 29-06-2015

%% Check if inputs exist
if nargin < 2 || isempty(isBenchmark)
    isBenchmark = 0; % input is real data
end

if nargin < 1 || isempty(input)
    error('Error: Please input a matrix');
end

%% Converting the input into all the formats

if size(input, 1) == size(input, 2) % If matrix
    Mat = input;
    numnodes = size(Mat, 1); % Number of nodes
    
    DirList = Mat2Direct(Mat, numnodes); % Calls function to convert matrix to directed list
    Undir = Direct2Undir(DirList); % Calls function to make undirected input
    
elseif size(input, 3) % List format
    if sum(input(:, 1) > input(:, 2)) == 0 % If undirected
        Undir = input; 
        numnodes = max(max(Undir(:, 1:2))); % Calculates the number of nodes in the system
        
        DirList = Undir2Direct(Undir); % Converts undirected to directed list
        Mat = Direct2Matrix(DirList, numnodes); % Calls function to make matrix input
    else % If directed
        
        DirList = input; 
        numnodes = max(max(DirList(:, 1:2))); % Calculates the number of nodes in the system
        
        Mat = Direct2Matrix(DirList, numnodes); % Calls function to make matrix input
        Undir = Direct2Undir(DirList); % Calls function to make undirected input
    end
else
    error('Error: Input is not one of the accepted formats');
end

%% Creating the final structure
Final = struct('Date', date);
Final.Network = Mat; % Saves the matrix of the network

%% Benchmark
if isBenchmark
    % Processes the data from the benchmark
    [BenchComm] = process_Benchmark(numnodes);
    
    % Places it in the final structure data
    Final.Benchmark = struct('Name', 'Benchmark', 'Result', BenchComm);
end

%% Running Functions
Final.Clauset = call_Clauset(DirList, numnodes);

for prec = [0 0.1 0.2 0.3]
    Final.(sprintf('Gopalan_prec_%g', prec*100)) = ...
        call_Gopalan(Undir, numnodes, 1000, prec);
end

Final.Jerry = call_Jerry(Mat, numnodes, 120, 0.9, 1, 'majority');

for prec = [0 0.01 0.1]
    Final.(sprintf('Link_prec_%g', prec*100)) = ...
        call_Link(Undir, numnodes, prec);
end

for thresh = [0 0.01 0.1]
Final.(sprintf('NNMF_thresh_%g', thresh*100)) = call_NNMF(Mat, thresh);
end

for tol = 0.1:0.1:1
Final.(sprintf('OSLOM_thresh_%g', tol*100)) = ...
    call_OSLOM(Undir, numnodes, 100, tol);
end

Final.Shen = call_Shen(Mat, numnodes, 7, 1.3, 0);

%% Saving data
fileName = 'Computation_Result.mat';
save(fileName, 'Final');

fprintf('All computation is complete! Saved as %s\n',fileName);

end