function [Output] = call_infomap(dirlist, numnodes,sourceCodePath)
% Infomap
% Input is Edgelist
% Runs the Infomap method

run_infomap(dirlist,sourceCodePath);
% 
%Processes the output
[Infomap_final] = process_infomap(numnodes);
% Places the structure data in
    Output= struct('Name', 'Infomap', 'Result', Infomap_final);
