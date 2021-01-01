function [Output] = call_infomap(dirlist, numnodes)
% Jerry
% Input is Edgelist
% Runs the Jerry method

run_infomap(dirlist);
% 
%Processes the output
[Infomap_final] = process_infomap(numnodes);
% Places the structure data in
    Output= struct('Name', 'Infomap', 'Result', Infomap_final);
