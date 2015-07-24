function Final = Computation(input, Methods,isBenchmark, benchfilename)
% Program that computes algorithms for a specific network
%
% INPUTS:
% input: adjacency matrix (or list of links in format ''node1, node2, weight'' as matrix)
% Methods: cell of names of methods to be tested
% isBenchmark: 1 for benchmark (otherwise 0)
% benchfilename: name of text file containing benchmark data
% Brandon Lam, 29-06-2015

%% IMPORTANT
% Before running this program, there are a few libraries that need to be
% loaded

% Gopalan - gsl and gcc libraries need to be loaded
% Link - requires python to be loaded

%% Check if inputs exist
if nargin < 1 || isempty(input)
    error('Error: Please input a matrix');
end

if nargin < 2 || isempty(Methods)
    Methods = {'Jerry', 'Shen'};
end

if nargin < 3 || isempty(isBenchmark)
    isBenchmark = 0; % input is real data
end

if isBenchmark == 1 && isempty(benchfilename)
    error('Error: Please input filename of benchmark community file');
end

%-------------------------------------------------------------------------------
%% Converting the input into all the formats
%-------------------------------------------------------------------------------

if size(input, 1) == size(input, 2) % If matrix
    Mat = input;
    numnodes = size(Mat, 1); % Number of nodes

    DirList = Mat2Direct(Mat, numnodes); % Calls function to convert matrix to directed list
    Undir = Mat2Undir(Mat); % Calls function to make undirected input

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
        Undir = Mat2Undir(Mat); % Calls function to make undirected input
    end
else
    error('Error: Input is not one of the accepted formats');
end

% So now we have 3 representations of the same object:
 % Mat: adjacency matrix (full)
 % DirList: list of edges
 % Undir: **NOT YET** undirected transformation of matrix (either edge exists counted as link)

%-------------------------------------------------------------------------------
%% Creating the final structure
% -- all the results of the specified algorithms is stored in a single structure
% called 'Final' (along with benchmark if specified, and data of run), and adjacency matrix
%-------------------------------------------------------------------------------
Final = struct('Date', date);
Final.Network = Mat; % Saves the matrix of the network

%% Benchmark
if isBenchmark
    % Processes the data from the benchmark
    [BenchComm] = process_Benchmark(benchfilename, numnodes);

    % Places it in the final structure data
    Final.Benchmark = struct('Name', 'Benchmark', 'Result', BenchComm);
end

%% Running Functions
for name = Methods
    switch name{1}
        case 'Clauset'
            Final.Clauset = call_Clauset(DirList, numnodes);
        case 'Gopalan'
            for prec = [0 0.1 0.2 0.3]
                Final.(sprintf('Gopalan_prec_%g', prec*100)) = ...
                    call_Gopalan(Undir, numnodes, 1000, prec);
            end
        case 'Jerry'
            Final.Jerry = call_Jerry(Mat, numnodes, 120, 0.09, 1, 'probabilistic');
        case 'Link'
            for prec = [0 0.01 0.1]
                Final.(sprintf('Link_prec_%g', prec*100)) = ...
                    call_Link(Undir, numnodes, prec);
            end
        case 'NNMF'
            for thresh = [0 0.01 0.1]
                Final.(sprintf('NNMF_thresh_%g', thresh*100)) = ...
                    call_NNMF(Mat, thresh);
            end
        case 'OSLOM'
            for tol = 0.1:0.1:1
                Final.(sprintf('OSLOM_thresh_%g', tol*100)) = ...
                    call_OSLOM(Undir, numnodes, 100, tol);
            end
        case 'Shen'
            Final.Shen = call_Shen(Mat, numnodes, 7, 1.3, 0);
    end
end

%% Saving data
fileName = 'Computation_Result.mat';
save(fileName, 'Final');

fprintf('All computation is complete! Saved as %s\n',fileName);

end
