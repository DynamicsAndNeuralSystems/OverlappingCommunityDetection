function[CommMat] = Clauset_rect_converter(~, ~)
% Converts the output of the Clauset method into a form the visualiser can use
% Brandon Lam, 23-07-2014

% Opens file, reads it.
fid = fopen('Clauset.txt');
Clausetdata = textscan(fid, '%f');
fclose(fid); % This imports the dat file into MATLAB
ClausetLabel = cell(length(Clausetdata{1}),1);
for Clau = 1:length(Clausetdata{1})
    ClausetLabel{Clau} = Clausetdata{1}(Clau);
end

% Replaces the output with communities that are ordered in size
ClausetLabel = comm_replace(ClausetLabel);

% Turns that into a matrix
CommMat{1} = NodeLabel2Mat(ClausetLabel);