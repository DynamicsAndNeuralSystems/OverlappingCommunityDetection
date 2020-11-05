function NMI = NMI_calc(Input1, Input2)
% Calculates the normalized similarity of two community assignments.
% **[NON-OVERLAPPING]**
%-------------------------------------------------------------------------------

numNodes = size(Input1, 1); % Calculates the number of nodes

% Binarize:
Comm1 = (Input1 ~= 0);
Comm2 = (Input2 ~= 0);

N = zeros(size(Comm1,2),size(Comm2,2));

for i = 1:numNodes
    % Creates the confusion matrix, where each row is a community of the
    % first input, and every column is the community of the second input.
    % each value is the number of nodes the input classifies as in that
    % specific community
    N(Comm1(i, :), Comm2(i, :)) = N(Comm1(i, :), Comm2(i, :)) + 1;
end

IXY = 0; % Initialising the for loops

% Calculating the IXY value
for x = 1:size(N,1)
    for y = 1:size(N,2)

        % P(x) = n(x)/n - First input
        Px = sum(N(x,:))/numNodes;
        % P(y) = n(y)/n - Second input
        Py = sum(N(:,y))/numNodes;
        % P(x,y) = n(x,y)/n - Overlap
        Pxy = N(x,y)/numNodes;

        if Pxy ~= 0
            % I(X,Y) = P(x,y)*log(P(x,y)/(P(x)*P(y)))
            IXY = IXY + Pxy*log(Pxy/(Px*Py));
        end
    end
end

HX = 0; % Initializing the next for loop

% Calculating H(X)
for x = 1:size(N, 1)

    % P(x) = n(x)/n - First input
    Px = sum(N(x,:))/numNodes;

    if Px ~= 0
        % H(X) = -P(x)*log(P(x))
        HX = HX - Px*log(Px);
    end
end

HY = 0; % Initialising the next for loop

% Calculating the H(Y) value
for y = 1:size(N, 2)

    % P(y) = n(y)/n - First input
    Py = sum(N(:, y))/numNodes;

    if Py ~= 0
        % H(Y) = -P(y)*log(P(y))
        HY = HY - Py*log(Py);
    end
end

% Final calculation
if (HX + HY) ~= 0
    NMI = 2*IXY/(HX + HY);
else
    NMI = 0;
end

end
