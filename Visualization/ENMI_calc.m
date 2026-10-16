function ENMI = ENMI_calc(Input1, Input2)
% Calculates the extended normalized value to calculate how similar two
% overlapping community matrices are.
% 
% Equations taken from Lancichinetti (2009). Detecting the overlapping and
% hierarchical community structure in complex networks.
%
% Inputs are numNodes x numComms membership matrices (nonzero = member; NaN
% is treated as 0).
% Note (a property of the LFK definition): a community containing all nodes or
% no nodes has zero entropy, and is scored as unmatched.
%-------------------------------------------------------------------------------

numNodes = size(Input1, 1); % Calculates the number of nodes

% Convert to binary:
Comm1 = (Input1 ~= 0) & ~isnan(Input1);
Comm2 = (Input2 ~= 0) & ~isnan(Input2);

%% Finding H(boldX|boldY), normalized
[HXXYYnorm] = extendNMI_calcs(Comm1, Comm2, numNodes);

%% Finding H(boldY|boldX), normalized
[HYYXXnorm] = extendNMI_calcs(Comm2, Comm1, numNodes);

%% Final calculation
ENMI = 1 - (HXXYYnorm + HYYXXnorm)/2;

end
