function Hnormalised = extendNMI_calcs(Comm1, Comm2, numNodes)
% This function calculates the H(boldX|boldY) normalised from the two
% inputs. This constitutes just under half of the calculations of the
% extended NMI.
% lower Hnormalised means better fit
% Brandon Lam, 16-07-2015
%-------------------------------------------------------------------------------

% H(boldX|boldY)
HXXYY = 0; % Initialising for loop

for k = 1:size(Comm1, 2)

    % H(Xk | Y) < Bold Y
    HXkYY = 10000; % Number to initiate the for loop, we want the smallest one possible
    for l = 1:size(Comm2, 2)

        PXk1Yl1 = sum((Comm1(:, k) + Comm2(:, l)) == 2)/numNodes; % P(Xk = 1, Yl = 1)
        PXk1Yl0 = sum(Comm1(:, k))/numNodes - PXk1Yl1; % P(Xk = 1, Yl = 0)
        PXk0Yl1 = sum(Comm2(:, l))/numNodes - PXk1Yl1; % P(Xk = 0, Yl = 1)
        PXk0Yl0 = 1 - (PXk0Yl1 + PXk1Yl0 + PXk1Yl1); % P(Xk = 0, Yl = 0)

        PYl1 = sum(Comm2(:, l))/numNodes; % P(Yl = 1)
        PYl0 = 1 - PYl1; % P(Yl = 0)

        % Constraint values - need them to COMPare against each other
        Comp1 = - PXk1Yl1*log(PXk1Yl1) - PXk0Yl0*log(PXk0Yl0);
        Comp2 = - PXk0Yl1*log(PXk0Yl1) - PXk1Yl0*log(PXk1Yl0);

        % Makes sure none of them are NaN, since that is a problem with
        % logs of 0
        if isnan(Comp1); Comp1 = 0; end
        if isnan(Comp2); Comp2 = 0; end

        if  Comp1 > Comp2
            % Constraint to protect from negative numbers

            % Temporary calculations, due to NaNs
            temp1 = - PXk1Yl1*log(PXk1Yl1);
            temp2 = - PXk0Yl0*log(PXk0Yl0);
            temp3 = - PXk0Yl1*log(PXk0Yl1);
            temp4 = - PXk1Yl0*log(PXk1Yl0);

            %% NaN fixing for HXkYl calculations

            if isnan(temp1); temp1 = 0; end
            if isnan(temp2); temp2 = 0; end
            if isnan(temp3); temp3 = 0; end
            if isnan(temp4); temp4 = 0; end

            %% Temp calcs for HYl, and NaN fixing

            temp5 = - PYl1*log(PYl1);
            temp6 = - PYl0*log(PYl0);
            if isnan(temp5); temp5 = 0; end
            if isnan(temp6); temp6 = 0; end

            %%

            HXkYl = temp1 + temp2 + temp3 + temp4; % H(Xk, Yl)
            HYl =  temp5 + temp6; % H(Yl)

            HXkgivenYl = HXkYl - HYl; % H(Xk|Yl)

            % Replaces the minimum number, to find the smallest one within this
            % loop
            if HXkgivenYl < HXkYY
                HXkYY = HXkgivenYl; % H(boldX|boldY)
            end
        end
    end

    PXk1 = sum(Comm1(:, k))/numNodes; % P(Xk = 1)
    PXk0 = 1 - PXk1; % P(Xk = 0)

    %% HXk calc and NaN fixing
    temp1 = - PXk1*log(PXk1);
    temp2 = - PXk0*log(PXk0);
    if isnan(temp1); temp1 = 0; end
    if isnan(temp2); temp2 = 0; end

    %%
    HXk = temp1 + temp2; % H(Xk)

    if HXkYY == 10000 % If no values fit the constraint
        HXkYYnorm = 1; % H(Xk|boldY) = H(Xk)
    else
        HXkYYnorm = HXkYY/HXk; % H(Xk|boldY) normalised
    end

    HXXYY = HXXYY + HXkYYnorm; % Adding to the total
end

Hnormalised = HXXYY/size(Comm1, 2); % H(boldX|boldY) normalised
