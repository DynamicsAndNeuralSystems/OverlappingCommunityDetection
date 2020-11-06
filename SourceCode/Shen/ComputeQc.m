function Qc = ComputeQc(Cover,adjMat)
% Computes the Q_c value for a given partition of nodes into groups
% Definition is Eq. (2) of Shen et al., 2009
%-------------------------------------------------------------------------------
% Ben Fulcher, 2014-03-31
%-------------------------------------------------------------------------------

if any(adjMat(adjMat>0)~=1)
    warning('The adjacency matrix is not binary... Tread carefully, padawan.');
end

if any(diag(adjMat)) ~= 0
    error('We need zeros on our diagonal.');
end

if sum(sum(adjMat-adjMat'))~=0
    error('Adjacency matrix should be symmetric');
end

fprintf(1,'This is a symmetric network with zero diagonals...! Well done.\n');

% adjMat is the adjacency matrix: [numNodes x numNodes]
% Cover is binary: [numNodes x numComms]

% ------------------------------------------------------------------------------
% Preliminaries
% ------------------------------------------------------------------------------
numNodes = size(Cover,1);
numComms = size(Cover,2);

% ------------------------------------------------------------------------------
% Compute belonging coefficients
% alpha: [numNodes x numComms]
% ------------------------------------------------------------------------------
alpha = ComputeBelongingCoeff(Cover,adjMat);

% ------------------------------------------------------------------------------
% Degree of each node
% k: [numNodes x 1]
% ------------------------------------------------------------------------------
k = sum(adjMat)';

% ------------------------------------------------------------------------------
% Compute Q_c
% ------------------------------------------------------------------------------

% Total weight of all edges, L
L = sum(k);

Qc_c = zeros(numComms,1);
for i = 1:numComms
    % Calculate the contribution to the sum from each community
    Qc_c(i) = sum(sum(alpha(:,i)*alpha(:,i)'.*(adjMat-k*k'/L)));
end
Qc = sum(Qc_c)/L;

fprintf(1,'Your Q_c value is %g, I reckon.\n',Qc);

end
