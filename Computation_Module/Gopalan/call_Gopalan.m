function [Output] = call_Gopalan(Undir, numnodes, maxIter, prec)
% Gopalan
% Input is either undirected or directed list, but it will treat them as
% undirected

% Check for inputs
if nargin < 3 || isempty(maxIter)
    maxIter = 1000; % Maximum iterations done within the algorithm
end
if nargin < 4 || isempty(prec)
    prec = 0.1; % Percentage of lowest links to cut before computation
end

% Gets rid of all links weaker than the percentage of precision
Input = Undir(Undir(:,3)>quantile(Undir(:,3),prec), :);

% Runs the Gopalan method
run_Gopalan(Input, numnodes, maxIter);

% Processes the data
[Gopalan_final] = process_Gopalan(numnodes);

% Puts the structural data in
Output = struct('Name', 'Gopalan', 'Percent_removed', prec*100, ...
    'Max_Iterations', maxIter, 'Result', Gopalan_final);
