function [Final,computationTimes] = OLCD_Compute(networkAdj,methodList,isBenchmark,benchmarkFileName,sourceCodePath,outputPath)
% Computes a set of (overlapping) community-detection algorithms for a given input network
%
%---INPUTS:
% networkAdj: adjacency matrix, or an edge list as a 3-column matrix
%             [node1, node2, weight] (node ids must be positive integers).
%             An edge list is treated as undirected if node1 <= node2 for all
%             rows, otherwise as directed.
%             NOTE: a 3x3 input is ambiguous (a 3-edge list or a 3-node matrix)
%             and is always treated as an adjacency matrix (with a warning).
% methodList: cell of strings labeling the methods to be evaluated (default
%             {'SLPA','Clique'}). Valid names (aliases in brackets):
%             'Infomap', 'OSLOM', 'NNMF', 'Link', 'SLPA' ('Jerry'),
%             'Clique' ('Shen'), 'WSBM' ('Clauset'), 'SVINET' ('Gopalan').
% isBenchmark: set true for a benchmark network (default false)
% benchmarkFileName: text file of benchmark node annotations (required if
%                    isBenchmark; each line is "node comm1 comm2 ...")
% sourceCodePath: path to the SourceCode directory (default: SourceCode in the
%                 toolbox root)
% outputPath: directory in which Computation_Result.mat is saved (default pwd)
%
%---OUTPUTS:
% Final: structure with the network, the benchmark (if any), and one field per
%        method run (e.g., Infomap, OSLOM_50, NNMF_10, Clique_5, Jerry, ...)
% computationTimes: structure with the computation time (seconds) of each
%                   method field in Final
%
%---IMPORTANT NOTE:
% Before running, some libraries must be available:
%   SVINET (Gopalan) - gsl and gcc libraries need to be loaded
%   Link - requires python to be loaded
%   OSLOM, Infomap - binaries must be built in SourceCode (see compile_all.sh)

%-------------------------------------------------------------------------------
%% Check Inputs
%-------------------------------------------------------------------------------
if nargin < 1 || isempty(networkAdj)
    error('Please input a matrix');
end
if nargin < 2 || isempty(methodList)
    methodList = {'Jerry','Shen'};
end
if ischar(methodList) || isstring(methodList)
    methodList = cellstr(methodList);
end
if nargin < 3 || isempty(isBenchmark)
    isBenchmark = false; % input is real data
end
if nargin < 4
    benchmarkFileName = '';
end
if isBenchmark && isempty(benchmarkFileName)
    error('Please input filename of benchmark community file');
end
if nargin < 5 || isempty(sourceCodePath)
    sourceCodePath = fullfile(fileparts(fileparts(mfilename('fullpath'))),'SourceCode');
end
if nargin < 6 || isempty(outputPath)
    outputPath = pwd;
end

% Validate method names (aliases map to canonical names):
aliasNames = {'Jerry','SLPA','Shen','Clique','Clauset','WSBM','Gopalan','SVINET','Infomap','OSLOM','NNMF','Link'};
canonicalNames = {'Jerry','Jerry','Shen','Shen','Clauset','Clauset','Gopalan','Gopalan','Infomap','OSLOM','NNMF','Link'};
[isValid,ix] = ismember(methodList,aliasNames);
if ~all(isValid)
    error('Unknown method(s): %s. Valid names: %s',strjoin(methodList(~isValid),', '),strjoin(aliasNames,', '));
end
methodList = canonicalNames(ix);

%-------------------------------------------------------------------------------
%% Converting the input into all the formats
%-------------------------------------------------------------------------------
if size(networkAdj,1) == size(networkAdj,2) % square adjacency matrix
    if size(networkAdj,1) == 3
        warning('A 3x3 input is ambiguous (3-edge list or 3-node matrix): treating as an adjacency matrix');
    end
    Mat = networkAdj;
    numNodes = size(Mat,1);

    DirList = Mat2Direct(Mat,numNodes); % directed edge list
    Undir = Mat2Undir(Mat); % undirected edge list

elseif size(networkAdj,2) == 3 % edge list
    nodeIds = networkAdj(:,1:2);
    if any(nodeIds(:) < 1) || any(nodeIds(:) ~= round(nodeIds(:)))
        error('Edge-list node ids must be positive integers');
    end
    numNodes = max(nodeIds(:));
    if all(networkAdj(:,1) <= networkAdj(:,2)) % undirected
        Undir = networkAdj;
        DirList = Undir2Direct(Undir);
        Mat = Direct2Matrix(DirList,numNodes);
    else % directed
        DirList = networkAdj;
        Mat = Direct2Matrix(DirList,numNodes);
        Undir = Mat2Undir(Mat);
    end
else
    error('Input network is not one of the accepted formats');
end

% So now we have 3 representations of the same network:
% Mat: adjacency matrix (full)
% DirList: directed edge list
% Undir: undirected edge list (weights (Aij+Aji)/2)

%-------------------------------------------------------------------------------
%% Creating the final structure
%-------------------------------------------------------------------------------
Final = struct('Date',date);
Final.Network = Mat;
computationTimes = struct();

%-------------------------------------------------------------------------------
%% Benchmark
%-------------------------------------------------------------------------------
if isBenchmark
    BenchComm = process_Benchmark(benchmarkFileName,numNodes);
    Final.Benchmark = struct('Name','Benchmark','Result',BenchComm);
end

%-------------------------------------------------------------------------------
%% Running Functions
%-------------------------------------------------------------------------------
numMethods = length(methodList);
fprintf(1,'Looping across %u OCDA methods.\n',numMethods);
for m = 1:numMethods
    switch methodList{m}
        case 'Clauset'
            runMethod('Clauset',@()call_Clauset(DirList,numNodes));
        case 'Infomap'
            runMethod('Infomap',@()call_infomap(DirList,numNodes,sourceCodePath));
        case 'Gopalan'
            for prec = [0,0.1,0.2,0.3]
                runMethod(sprintf('Gopalan_prec_%d',round(prec*100)), ...
                    @()call_Gopalan(Undir,numNodes,1000,prec,sourceCodePath));
            end
        case 'Jerry'
            runMethod('Jerry',@()call_Jerry(Mat,numNodes,120,0.09,1,'probabilistic'));
        case 'Link'
            for prec = [0,0.01,0.1]
                runMethod(sprintf('Link_prec_%d',round(prec*100)), ...
                    @()call_Link(Undir,numNodes,prec,sourceCodePath));
            end
        case 'NNMF'
            for thresh = [0.1,0.2,0.3,0.4]
                runMethod(sprintf('NNMF_%d',round(thresh*100)),@()call_NNMF(Mat,thresh));
            end
        case 'OSLOM'
            % (tol = 1 is degenerate: every module passes significance)
            for tol = 0.1:0.1:0.9
                runMethod(sprintf('OSLOM_%d',round(tol*100)), ...
                    @()call_OSLOM(Undir,numNodes,100,tol,sourceCodePath));
            end
        case 'Shen'
            for cliqueSize = [3,4,5,6,7,9]
                runMethod(sprintf('Clique_%d',cliqueSize), ...
                    @()call_Shen(0.5*(Mat+Mat'),numNodes,cliqueSize,1.3,0));
            end
    end
end

%-------------------------------------------------------------------------------
%% Saving data
%-------------------------------------------------------------------------------
fileName = fullfile(outputPath,'Computation_Result.mat');
save(fileName,'Final','computationTimes');
fprintf('All computation is complete! Saved as %s.\n',fileName);

%-------------------------------------------------------------------------------
    function runMethod(fieldName,callFunction)
        % Runs a method, storing its result and computation time
        methodTimer = tic;
        Final.(fieldName) = callFunction();
        computationTimes.(fieldName) = toc(methodTimer);
    end

end
