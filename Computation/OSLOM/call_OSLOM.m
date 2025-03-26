function Output = call_OSLOM(Undir, numNodes, numIter, Tol, sourceCodePath, outputPath)
%% OSLOM
% Input, Undir, is a sparse matrix (undirected)
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs
%-------------------------------------------------------------------------------
if nargin < 3 || isempty(numIter)
    numIter = 100; % Number of iterations within the algorithm
end
if nargin < 4 || isempty(Tol)
    Tol = 0.5; % Range of tolerances
end
%-------------------------------------------------------------------------------

% Run OSLOM on the matrix:
run_OSLOM(Undir, numIter, Tol, sourceCodePath, outputPath);

% Process the text file outputs:
OSLOM_final = process_OSLOM(Tol,numNodes, outputPath);

% Form an output structure:
Output = struct('Name','OSLOM','Threshold',Tol,'Result',OSLOM_final{1});

end
