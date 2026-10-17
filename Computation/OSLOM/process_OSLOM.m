function Output = process_OSLOM(tpFile,numNodes)
% Processes the OSLOM module file (one module per line, listing node ids) into
% a numNodes x numModules membership matrix.
%-------------------------------------------------------------------------------

modules = readCommunityLines(tpFile);
Output = communitiesToMatrix(modules,numNodes);

end
