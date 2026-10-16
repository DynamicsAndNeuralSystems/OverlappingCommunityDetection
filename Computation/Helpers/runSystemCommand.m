function cmdOut = runSystemCommand(cmd)
% Runs a shell command, raising an error (with the command's output) if it fails,
% rather than silently continuing and later reading stale or missing output.
%-------------------------------------------------------------------------------

[status,cmdOut] = system(cmd);
if status ~= 0
    error('Command failed (exit status %d):\n%s\n%s',status,cmd,cmdOut);
end

end
