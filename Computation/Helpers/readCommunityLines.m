function lineValues = readCommunityLines(fileName)
% Reads a text file of whitespace-separated integers, one community (or node)
% per line, skipping blank lines and comment lines beginning with '#'.
% Returns a cell array with one (column) vector of integers per line.
% Does not rely on trailing whitespace (unlike textscan-based parsing).
%-------------------------------------------------------------------------------

if ~isfile(fileName)
    error('Cannot find the output file %s',fileName);
end
textLines = splitlines(fileread(fileName));
textLines = strtrim(textLines);
textLines = textLines(~cellfun(@isempty,textLines) & ~startsWith(textLines,'#'));
lineValues = cellfun(@(x)sscanf(x,'%f'),textLines,'UniformOutput',false);

end
