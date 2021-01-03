function F1 = F1_overlapcalc(Input1, Input2)
% Calculates the F1 score between the true and predicted set of overlapping nodes.
% Input1: True community structure
% Input2: Community structure predicted by an OCDA
% 
% F1 = 2*(precision*recall)/(precision+recall)
%-------------------------------------------------------------------------------
% Extracting set of overlapping nodes

binInput1 = Input1>0;
% Number of communities assigned to a node
ncomm1 = sum(binInput1,2);
% Logical describing whether a node is overlapping or not
overlap1 = ncomm1>1;

binInput2 = Input2>0;
% Number of communities assigned to a node
ncomm2 = sum(binInput2,2);
% Logical describing whether a node is overlapping or not
overlap2 = ncomm2>1;

% Including this condition to avoid NaNs and to reflect the method's
% inability to identify overlapping nodes
if sum(overlap1)==0 | sum(overlap2)==0
    F1 = 0;
else
    % Calling function to compute binary F1 score
    F1 = F1_calc(overlap1, overlap2);
end