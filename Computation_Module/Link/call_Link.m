function [Output] = call_Link(Undir, numnodes, prec)
% Link
% Input is undirected list
% Check for inputs
if nargin < 3 || isempty(prec)
    prec = 0.01; % Percentage of lowest links to cut before computation
end

    % Gets rid of all links weaker than the percentage of precision
    Input = Undir(Undir(:,3)>quantile(Undir(:,3),prec), :);
    
    % Runs the link communities method
    run_Link(Input);
    
    % Processes the textfile outputs
    Link_final = process_Link(numnodes);
    
    % Puts the structural data in
    Output = ...
        struct('Name', 'Link', 'Percent_removed', prec*100, ...
        'Result', Link_final);