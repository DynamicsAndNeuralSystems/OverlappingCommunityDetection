function NMI = NMI_calc_vecs(Input1, Input2)
% Calculates the normalized similarity of two community assignments.
% **[NON-OVERLAPPING]**
% Inputs are vectors (row or column) of positive integer community labels.
%-------------------------------------------------------------------------------

Input1 = Input1(:);
Input2 = Input2(:);
numNodes = numel(Input1); % Calculates the number of nodes

% Get the number of unique communities
numComm1 = max(Input1);
numComm2 = max(Input2);

% Compute confusion matrix N
N = accumarray([Input1, Input2], 1, [numComm1, numComm2]);

% Compute joint probability P(X,Y)
Pxy = N / numNodes;
% Compute marginal probabilities P(X) and P(Y)
Px = sum(Pxy, 2); % Sum over columns
Py = sum(Pxy, 1); % Sum over rows

% Compute mutual information I(X, Y) (summing over nonzero joint probabilities)
PxPy = Px*Py;
isPos = (Pxy > 0);
IXY = sum(Pxy(isPos) .* log(Pxy(isPos) ./ PxPy(isPos)));

% Compute entropies H(X) and H(Y)
HX = -sum(Px(Px > 0) .* log(Px(Px > 0)));
HY = -sum(Py(Py > 0) .* log(Py(Py > 0)));

% Compute Normalized Mutual Information
if (HX + HY) > 0
    NMI = 2 * IXY / (HX + HY);
else
    NMI = 0;
end

end
