function ENMI = ENMI_calc(Input1, Input2)
% Calculates the extended normalized value to calculate how similar two
% overlapping community matrices are.
% 
% Equations taken from Lancichinetti (2009). Detecting the overlapping and
% hierarchical community structure in complex networks.
%-------------------------------------------------------------------------------

numNodes = size(Input1, 1); % Calculates the number of nodes

% Convert to binary:
Comm1 = (Input1 ~= 0);
Comm2 = (Input2 ~= 0);

%% Finding H(boldX|boldY), normalized
[HXXYYnorm] = extendNMI_calcs(Comm1, Comm2, numNodes);

%% Finding H(boldY|boldX), normalized
[HYYXXnorm] = extendNMI_calcs(Comm2, Comm1, numNodes);

%% Final calculation
ENMI = 1 - (HXXYYnorm + HYYXXnorm)/2;

end
