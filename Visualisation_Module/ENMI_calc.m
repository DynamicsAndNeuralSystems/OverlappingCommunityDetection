function [ENMI] = ENMI_calc(Input1, Input2)
% Calculates the extended normalised value to calculate how similar two
% overlapping community matrices are.
% Brandon Lam, 16-07-2015

numnodes = size(Input1, 1); % Calculates the number of nodes

Comm1 = Input1 ~= 0; % Binary Input
Comm2 = Input2 ~= 0; % Binary Input

% For all equations, look to Lancichinetti, 2009 - Detecting the
% overlapping and heirarchical community structure in complex networks

%% Finding H(boldX|boldY) normalised
[HXXYYnorm] = extendNMI_calcs(Comm1, Comm2, numnodes);

%% Finding H(boldY|boldX) normalised
[HYYXXnorm] = extendNMI_calcs(Comm2, Comm1, numnodes);

%% Final calculations
ENMI = 1 - (HXXYYnorm + HYYXXnorm)/2;
