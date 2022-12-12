function [ssMod] = modifiedLayup(ss)

% If even number of plies
if rem(length(ss), 2) == 0

    % Add 'extra' ply angle to middle
    % Extra ply is left of the midplane
    ssMod = [ss(1: length(ss)/2), ss(length(ss)/2), ss((length(ss)/2)+1: end)];
    
% If odd number of plies
else
    
    % Add 'extra' ply angle to middle
    % Extra ply is the middle ply
    ssMod = [ss(1:ceil(length(ss)/2)), ss(ceil(length(ss)/2)), ...
        ss(ceil(length(ss)/2) + 1: (length(ss)))];
end

end
