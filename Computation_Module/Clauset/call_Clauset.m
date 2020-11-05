function Output = call_Clauset(DirList, numnodes)
% Clauset - uses weighted directed list as an input
% This method is relatively easy - it uses a directed list of connections
% as an input, but the output has no overlaps.

Clauset = run_Clauset(DirList); % Runs the Clauset method
Clauset_final = process_Clauset(Clauset, numnodes); % Processes the output

% Putting it in the structure
Output = struct('Name', 'Clauset', 'Result', Clauset_final);

end
