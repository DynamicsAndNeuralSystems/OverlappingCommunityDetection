function Hnormalised = extendNMI_calcs(Comm1, Comm2, numNodes)
% This function calculates the H(boldX|boldY) normalised from the two
% inputs. This constitutes just under half of the calculations of the
% extended NMI (Lancichinetti, Fortunato & Kertesz, 2009).
% lower Hnormalised means better fit
% Comm1, Comm2: binary numNodes x numComms membership matrices (NaN is
% treated as 0).
% Note (a property of the LFK definition): a community containing all nodes or
% no nodes has zero entropy, and is scored as unmatched (normalised value 1).
% Brandon Lam, 16-07-2015 (vectorised version of the original loops)
%-------------------------------------------------------------------------------

% Treat NaN as no membership; double for matrix products
Comm1 = double(Comm1 ~= 0 & ~isnan(Comm1));
Comm2 = double(Comm2 ~= 0 & ~isnan(Comm2));
numComm1 = size(Comm1, 2);
numComm2 = size(Comm2, 2);

% Pairwise probabilities, rows index X_k, columns index Y_l
PXk1 = (sum(Comm1, 1)'/numNodes) * ones(1, numComm2); % P(Xk = 1)
PYl1 = ones(numComm1, 1) * (sum(Comm2, 1)/numNodes); % P(Yl = 1)
PYl0 = 1 - PYl1; % P(Yl = 0)
PXk1Yl1 = (Comm1'*Comm2)/numNodes; % P(Xk = 1, Yl = 1)
PXk1Yl0 = PXk1 - PXk1Yl1; % P(Xk = 1, Yl = 0)
PXk0Yl1 = PYl1 - PXk1Yl1; % P(Xk = 0, Yl = 1)
PXk0Yl0 = 1 - (PXk0Yl1 + PXk1Yl0 + PXk1Yl1); % P(Xk = 0, Yl = 0)

% Entropy terms, with 0*log(0) = NaN replaced by 0
temp1 = nanTo0(-PXk1Yl1.*log(PXk1Yl1));
temp2 = nanTo0(-PXk0Yl0.*log(PXk0Yl0));
temp3 = nanTo0(-PXk0Yl1.*log(PXk0Yl1));
temp4 = nanTo0(-PXk1Yl0.*log(PXk1Yl0));
temp5 = nanTo0(-PYl1.*log(PYl1));
temp6 = nanTo0(-PYl0.*log(PYl0));

% Constraint values - need them to compare against each other
Comp1 = nanTo0(-PXk1Yl1.*log(PXk1Yl1) - PXk0Yl0.*log(PXk0Yl0));
Comp2 = nanTo0(-PXk0Yl1.*log(PXk0Yl1) - PXk1Yl0.*log(PXk1Yl0));

HXkYl = temp1 + temp2 + temp3 + temp4; % H(Xk, Yl)
HYl = temp5 + temp6; % H(Yl)
HXkgivenYl = HXkYl - HYl; % H(Xk|Yl)

% Only those pairs satisfying the constraint are candidates:
HXkgivenYl(~(Comp1 > Comp2)) = Inf;
HXkgivenYl(~(HXkgivenYl < 10000)) = Inf;
if numComm2 > 0
    HXkYY = min(HXkgivenYl, [], 2); % H(Xk|boldY)
else
    HXkYY = Inf(numComm1, 1);
end

% H(Xk)
PXk1 = sum(Comm1, 1)'/numNodes;
HXk = nanTo0(-PXk1.*log(PXk1)) + nanTo0(-(1 - PXk1).*log(1 - PXk1));

HXkYYnorm = HXkYY./HXk; % H(Xk|boldY) normalised
HXkYYnorm(isinf(HXkYY)) = 1; % If no values fit the constraint, H(Xk|boldY) = H(Xk)

Hnormalised = sum(HXkYYnorm)/numComm1; % H(boldX|boldY) normalised

end

function x = nanTo0(x)
% Makes sure none of them are NaN, since that is a problem with logs of 0
x(isnan(x)) = 0;
end
