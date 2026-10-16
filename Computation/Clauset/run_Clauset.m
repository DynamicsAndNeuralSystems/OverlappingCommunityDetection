function [Labels,Model] = run_Clauset(DirList,k)
% Function that runs the Clauset method (WSBM) on a directed list of connections
% with k blocks (default: 4, the wsbm default)

if nargin < 2 || isempty(k)
    k = 4;
end

% Runs the function of the Clauset algorithm (silently)
[Labels,Model] = wsbm(DirList,k,'verbosity',0);

end
