function [tempDir,cleanupObj] = makeTempDir()
% Creates a fresh temporary directory for an algorithm's input/output files, so
% that runs don't write into SourceCode/ and concurrent runs can't clobber each
% other's files. The directory is deleted when cleanupObj is cleared (e.g., when
% the calling function returns or errors).
%-------------------------------------------------------------------------------

tempDir = tempname;
mkdir(tempDir);
cleanupObj = onCleanup(@()rmdir(tempDir,'s'));

end
