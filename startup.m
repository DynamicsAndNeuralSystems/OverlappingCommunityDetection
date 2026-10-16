% Add paths required for the project (ignoring hidden folders, including version control).
% Uses the location of this file, so it works from any current directory.
rootDir = fileparts(mfilename('fullpath'));
files = dir(rootDir);
isDir = [files.isdir] & ~startsWith({files.name},'.'); % directories, removing hidden (and '.', '..')
dirNames = {files(isDir).name};
for j = 1:length(dirNames)
    pathList = strsplit(genpath(fullfile(rootDir,dirNames{j})),pathsep);
    pathList = pathList(~cellfun(@isempty,pathList));
    % Remove hidden subdirectories (e.g., .git):
    isHidden = contains(pathList,[filesep '.']);
    addpath(pathList{~isHidden})
end
