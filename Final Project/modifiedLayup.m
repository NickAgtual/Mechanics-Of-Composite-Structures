function [ssMod] = modifiedLayup(ss)

if rem(length(ss), 2) == 0
    
    ssMod = [ss(1: length(ss)/2), ss(length(ss)/2), ss(length(ss)+1: end)];
    
else
    
    x = 1;
end

end
