function Qc = ComputeQc(Cover,AdjMat)
% Computes the Q_c value for a given partition of nodes into groups
% Ben Fulcher, 2014-03-31

% Definition is Eq. (2) of Shen et al., 2009

if any(AdjMat(AdjMat>0)~=1)
    warning('The adjacency matrix is not binary... Tread carefully, padawan.');
end

if any(diag(AdjMat)) ~= 0
    error('We need zeros on our diagonal.');
end

if sum(sum(AdjMat-AdjMat'))~=0
    error('Adjacency matrix should be symmetric');
end

fprintf(1,'This is a symmetric network with zero diagonals...! Well done.\n');

% AdjMat is the adjacency matrix: [NumNodes x NumNodes]
% Cover is binary: [NumNodes x NumComms]

% ------------------------------------------------------------------------------
% Preliminaries
% ------------------------------------------------------------------------------
NumNodes = size(Cover,1);
NumComms = size(Cover,2);

% ------------------------------------------------------------------------------
% Compute belonging coefficients
% ------------------------------------------------------------------------------
alpha = ComputeBelongingCoeff(Cover,AdjMat);

% alpha: [NumNodes x NumComms]

% ------------------------------------------------------------------------------
% Degree of each node
% ------------------------------------------------------------------------------

k = sum(AdjMat)';

% k: [NumNodes x 1]

% ------------------------------------------------------------------------------
% Compute Q_c
% ------------------------------------------------------------------------------

% Total weight of all edges, L
L = sum(k);

Qc_c = zeros(NumComms,1);
for i = 1:NumComms
    % Calculate the contribution to the sum from each community
    Qc_c(i) = sum(sum(alpha(:,i)*alpha(:,i)'.*(AdjMat-k*k'/L)));
end
Qc = sum(Qc_c)/L;

fprintf(1,'Your Q_c value is %g, I reckon.\n',Qc);

end