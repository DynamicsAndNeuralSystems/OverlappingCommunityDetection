function NMI = NMI_calc_vecs(Input1, Input2)
% Calculates the normalized similarity of two community assignments.
% **[NON-OVERLAPPING]**
%-------------------------------------------------------------------------------

numNodes = size(Input1, 1); % Calculates the number of nodes

% Get the number of unique communities
numComm1 = max(Input1);
numComm2 = max(Input2);

% Compute confusion matrix N
N = accumarray([Input1(:), Input2(:)], 1, [numComm1, numComm2]);

% Compute joint probability P(X,Y)
Pxy = N / numNodes;
% Compute marginal probabilities P(X) and P(Y)
Px = sum(Pxy, 2); % Sum over columns
Py = sum(Pxy, 1); % Sum over rows

% Compute mutual information I(X, Y)
IXY = 0;
for x = 1:numComm1
    for y = 1:numComm2
        if Pxy(x, y) > 0
            IXY = IXY + Pxy(x, y) * log(Pxy(x, y) / (Px(x) * Py(y)));
        end
    end
end

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