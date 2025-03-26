function Final = OLCD_Compute(networkAdj,methodList,isBenchmark,benchmarkFileName,sourceCodePath, outputPath)
% Computes a set of (overlapping) community-detection algorithms for a given input network
%
%---INPUTS:
% input: adjacency matrix (or list of links in format ''node1, node2, weight'' as matrix)
% methodList: cell of strings labeling the methods to be evaluated
% isBenchmark: set true for a benchmark network
% benchmarkFileName: name of a text file containing benchmark node annotations
%
%%---IMPORTANT NOTE:
% Before running this program, there are a few libraries that need to be
% loaded:
%   Gopalan - gsl and gcc libraries need to be loaded
%   Link - requires python to be loaded

%-------------------------------------------------------------------------------
%% Check Inputs
%-------------------------------------------------------------------------------
if nargin < 1 || isempty(networkAdj)
    error('Error: Please input a matrix');
end

if nargin < 2 || isempty(methodList)
    methodList = {'Jerry', 'Shen'};
end

if nargin < 3 || isempty(isBenchmark)
    isBenchmark = false; % input is real data
end

if isBenchmark == 1 && isempty(benchmarkFileName)
    error('Error: Please input filename of benchmark community file');
end

%-------------------------------------------------------------------------------
%% Converting the input into all the formats
%-------------------------------------------------------------------------------

if size(networkAdj,1) == size(networkAdj,2) % square adjacency matrix
    Mat = networkAdj;
    numNodes = size(Mat,1); % Number of nodes

    DirList = Mat2Direct(Mat,numNodes); % Calls function to convert matrix to directed list
    Undir = Mat2Undir(Mat); % Calls function to make undirected networkAdj

elseif size(networkAdj,2)==3 % List format
    if sum(networkAdj(:,1) > networkAdj(:,2)) == 0 % (undirected)
        Undir = networkAdj;
        numNodes = max(max(Undir(:,1:2))); % Calculates the number of nodes in the system

        DirList = Undir2Direct(Undir); % Converts undirected to directed list
        Mat = Direct2Matrix(DirList,numNodes); % Calls function to make matrix networkAdj
    else % If directed
        DirList = networkAdj;
        numNodes = max(max(DirList(:,1:2))); % Calculates the number of nodes in the system

        Mat = Direct2Matrix(DirList,numNodes); % Calls function to make matrix input
        Undir = Mat2Undir(Mat); % Calls function to make undirected input
    end
else
    error('Error: Input adjacency matrix is not one of the accepted formats');
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

%-------------------------------------------------------------------------------
%% Benchmark
%-------------------------------------------------------------------------------
if isBenchmark
    % Processes the data from the benchmark
    BenchComm = process_Benchmark(benchmarkFileName,numNodes);

    % Places it in the final structure data
    Final.Benchmark = struct('Name','Benchmark','Result',BenchComm);
end

%-------------------------------------------------------------------------------
%% Running Functions
%-------------------------------------------------------------------------------


numMethods = length(methodList);
fprintf(1,'Looping across %u OCDA methods.\n',numMethods);
for m = 1:numMethods
    theMethod = methodList{m};
    switch theMethod
        case 'Clauset'
            Final.Clauset = call_Clauset(DirList, numNodes);
        case 'Infomap'
            Final.Infomap = call_infomap(DirList, numNodes, sourceCodePath);
        case 'Gopalan'
            for prec = [0 0.1 0.2 0.3]
                Final.(sprintf('Gopalan_prec_%g', prec*100)) = ...
                    call_Gopalan(Undir, numNodes, 1000, prec);
            end
        case 'Jerry'
            Final.Jerry = call_Jerry(Mat, numNodes, 120, 0.09, 1, 'probabilistic');
        case 'Link'
            for prec = [0,0.01,0.1]
                Final.(sprintf('Link_prec_%g', prec*100)) = ...
                    call_Link(Undir, numNodes, prec);
            end
        case 'NNMF'
            for thresh = [0.1,0.2,0.3,0.4]
                Final.(sprintf('NNMF_%g', thresh*100)) = ...
                    call_NNMF(Mat, thresh);
            end
        case 'OSLOM'
            for tol = 0.1:0.1:1
                Final.(sprintf('OSLOM_%g', tol*100)) = ...
                    call_OSLOM(Undir, numNodes, 100, tol, sourceCodePath, outputPath);
            end
        case 'Shen'
            for clique_size = [3,4,5,6,7,9]
                Final.(sprintf('Clique_%g', clique_size)) = ...
                    call_Shen(0.5*(Mat+Mat'), numNodes, clique_size, 1.3, 0);
            end
    end
end

%-------------------------------------------------------------------------------
%% Saving data
%-------------------------------------------------------------------------------
fileName = fullfile(outputPath, 'Computation_Result.mat');
save(fileName, 'Final');
fprintf('All computation is complete! Saved as %s.\n',fileName);

end
