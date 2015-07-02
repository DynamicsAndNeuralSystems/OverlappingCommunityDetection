function[Label] = comm_replace(Label)
% Turns a cell array of communities into the same thing, but with the
% communities ordered in size
% Brandon Lam, 14-07-2014

% Obtains the sorted communities in relation to each other
    [~, Replace] = sort(tabulate(horzcat(Label{:})));
    % Inverts it so we can convert the original easier
    [~, Inv] = sort(Replace(:, 3));
    % Replaces all community values with the new sorted ones
    for Ind = 1:size(Label, 1)
        % Works out which one it needs to replace it with, and sorts it so
        % the colours will appear on the same side, usually
        Label{Ind} = sort(Inv(Label{Ind})'); 
    end
end