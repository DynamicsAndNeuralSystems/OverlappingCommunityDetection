function [NodeLabels, NumComms] = Converter(file)

% Brandon Lam, 04-07-2014
% Converts the output file of the OSLOM function to a form that the brain
% visualiser takes.

fid = fopen(file); % Opens the inputted file
% Reads the text file, outputting all the numbers (where new lines
% separate them with a 0
A = textscan(fid, '%d', 'CommentStyle', '#', 'Delimiter', {' ', '\n'});
fclose(fid);

%% [~, Col] = size(A{1}'); % specifies the number of columns needed

NumComms = sum(A{1}==0); % Finds the total number of communities
NumNodes = max(A{1}); % Finds the number of nodes
NodeLabels = cell(NumNodes,1); % Using NumNodes, it makes cells for NodeLabels
Delims = [0, find(A{1}==0)']; % Gives the distances between each module change
% This for loop reads every single node's community, adding the module to
% the specific row's cell.
for i = 1:NumComms % For the number of communities
    for Pos = Delims(i)+1:Delims(i+1)-1 % For the nodes between the Delims
        Node = A{1}(Pos); % Finds the specific node number
        NodeLabels{Node} = [NodeLabels{Node}, i]; % Adds it to the cell
    end
end

%% Module = 1; % Allocating the first module
% Counter = 1; % Used to place numbers where they should be
% 
% % The purpose of this for loop is to convert the read file into a nice
% % matrix, where the rows are the modules, and the numbers inside are the
% % node IDs.
% for Check = 1:Col
%     if A{1}(Check, 1) == 0 % "If there is a 0, change to the next module"
%         Module = Module + 1;
%         Counter = 1; % Resets the counter to 1 so the next line starts on the first matrix place.
%     else
%         mtx(Module, Counter) = A{1}(Check, 1); % The number is then placed in the right spot in the matrix
%         Counter = Counter + 1; % Moves the counter over so the number is saved
%     end
% end
% 
% [Row, Col] = size(mtx); % Reallocating variables
% Arb = 1; % Arbitrary number, used to be the counter
% ColMat = zeros(max(A{1}), 1); % First reallocation of the new matrix
% 
% 
% % The purpose of these loops is to convert the above matrix into another
% % matrix, where the rows are the node IDs, and the values and the modules
% % it belongs to.
% for Module = 1:Row % For all numbers in the mtx matrix
%     for Nodes = 1:Col
%         if mtx(Module, Nodes) == 0
%             break % If there's a 0 in the column, then no more numbers are 
%                   % after this column. This break will cause it to go to 
%                   % the next module (more efficient)
%         else
%             while ColMat(mtx(Module, Nodes), Arb) ~= 0 % If the value is not 0, we need to move it over
%                 Arb = Arb + 1; % Add one to the arbitrary number to shift our focus one right
%                 [~, x] = size(ColMat); % Define new variable, allowing to check whether this new
%                                            % column has been added already
%                 if x < Arb % If it hasn't, add another one
%                     ColMat(max(A{1}), end+1) = 0; % Adds another column to the matrix
%                 end
%             end
%             ColMat(mtx(Module, Nodes), Arb) = Module; % Sets the value into the new matrix
%             Arb = 1; % Resets the arbitrary number back to 1, since we now have a different number
%         end
%     end
% end
% 
% NodeLabels = num2cell(ColMat, 2); % Output is the cells
% NumComms = Row; % Number of communities
% end