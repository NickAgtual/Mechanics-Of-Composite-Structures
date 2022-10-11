clear; clc; close all
%% Global System

globSys.e1 = [1 0 0]';
globSys.e2 = [0 1 0]';
globSys.e3 = [0 0 1]';

globSys.stress = [4 3 0 0 0 2]'; % MPa

%% Local System

locSys.e1 = [(sqrt(2)/2) (-sqrt(2)/2) 0]';
locSys.e2 = [(sqrt(2)/2) (sqrt(2)/2) 0]';
locSys.e3 = [0 0 1]';

%% Transformation Matrix

% Combining all unit vectors into a matrix
globSys.comb = [globSys.e1 globSys.e2 globSys.e3];
locSys.comb = [locSys.e1 locSys.e2 locSys.e3];

% Determining transformation matrix
for ii = 1: length(globSys.e1)
    for jj = 1: length(locSys.e1)
        
        b(ii, jj) = locSys.comb(jj, ii);
          
    end
end


%% Stress Transformation Matrix

stressTransform = [b(1, 1) ^ 2, b(1, 2) ^ 2, ...
    b(1, 3) ^ 2, 2 * b(1, 2) * b(1, 3), 2 * b(1, 1) * ...
    b(1, 3), 2 * b(1, 1) * b(1, 2); b(2, 1) ^ 2 ,...
    b(2, 2) ^ 2, b(2, 3) ^ 2, 2 * b(2, 2) * b(2, 3), 2 * b(2, 1) * ...
    b(2, 3), 2 * b(2, 1) * b(2, 2); b(3, 1) ^ 2, b(3, 2) ^ 2, ...
    b(3, 3) ^ 2, 2 * b(3, 2) * b(3, 3), 2 * b(3, 1) * b(3, 3), ...
    2 * b(3, 1) * b(3, 2); b(2, 1) * b(3, 1), b(2, 2) * b(3, 3), ...
    b(2, 3) * b(3, 3), (b(2, 2) * b(3, 3)) + (b(2, 3) * b(3, 2)), ...
    (b(2, 3) * b(3, 1)) + (b(2, 1) * b(3, 3)), (b(3, 1) * b(2, 2)) + ...
    (b(2, 1) * b(3, 2)); b(1, 1) * b(3, 1), b(1, 2) * b(3, 2), b(1, 3) ...
    * b(3, 3), (b(1, 3) * b(3, 2)) + (b(1, 2) * b(3, 3)), ...
    (b(1, 1) * b(3, 3)) + (b(1, 3) * b(3, 1)), (b(1, 1) * b(3, 2)) + ...
    (b(1, 2) * b(3, 1)); b(1, 1) * b(2, 1), b(1, 2) * b(2, 2), b(1, 3) ...
    * b(2, 3), (b(1, 2) * b(2, 3)) * (b(1, 3) * b(2, 2)), (b(1, 1) * ...
    b(2, 3)) + (b(1, 3) * b(2, 1)), (b(1, 1) * b(2, 2)) + (b(1, 2) * ...
    b(2, 1))];

locSys.stress = stressTransform * globSys.stress;
locSys.stressTensor = [locSys.stress(1) locSys.stress(4) locSys.stress(6); ...
    locSys.stress(4) locSys.stress(2) locSys.stress(5); ...
    locSys.stress(6) locSys.stress(5) locSys.stress(3)];

%% Stress Vector and Normal Components on Surface ABCD

% Each row is a stress vector
locSys.stressVec = zeros(size(locSys.comb, 2), size(locSys.comb, 1));

for ii = 1:size(locSys.comb, 2)
    locSys.stressVec(ii, :) = locSys.comb(:, ii)' * locSys.stressTensor;
    
    % Each element corresponds to each stress vec
    locSys.norm(ii) = locSys.stressVec(ii, :) * locSys.comb(ii,:)';
end







