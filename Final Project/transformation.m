function [locStress, globStress, locStrain, globStrain] = ...
    transformation(stressOrStrain, state, globOrLoc,  theta, moduli)

%% Function Definition

% Function Outputs:
% 1) Local state of stress
% 2) Global state of stress
% 3) Local state of strain
% 4) Global state of strain

% Function Inputs:
% 1) If input is state of stress, input 1, elseif input is state of strain
% input 0
% 2) State of stress or strain (3 element row vector or column vector)
% 3) If input is state is global, input 1, elseif input is state is local
% input 0
% 4) Angle between global and local coordinate system (+ CCW)

%% Specifying Input Parameters

% Checking whether state of stress/strain was inputted as row or col vec
if isequal(size(state), [1, 3])
    % If row vector, convert to column vector
    state = state';
elseif ~isequal(size(state), [3, 1], [1, 3])
    % Checking for invalid stress/strain input
    fprintf('Invalid State of Stress or Strain')
end

% Assigning user input to local/global stress/strain based on input
if stressOrStrain == 1 && globOrLoc == 1
    globSys.stress = state;
    condition = 1;
elseif stressOrStrain == 1 && globOrLoc == 0
    locSys.stress = state;
    condition = 2;
elseif stressOrStrain == 0 && globOrLoc == 1
    globSys.strain = state;
    condition = 3;
elseif stressOrStrain == 0 && globOrLoc == 0
    locSys.strain = state;
    condition = 4;
else
    fprintf('Invalid Input Arguments. See Function Definition for IO')
end

%% Reduced Compliance and Stiffness Matrices

% Reduced Compliance Matrix
S = [(1/moduli.E1) (-moduli.v12/moduli.E1) 0; ...
    (-moduli.v12/moduli.E1) (1/moduli.E2) 0; ...
    0 0 (1/moduli.G12)];

% Reduced stiffness matrix
Q = inv(S);

%% Transformed Reduced Stiffness and Compliance Matrices

% Reduced stress transformation
reducedStressTransformation = @(theta) ...
    [cosd(theta)^2 sind(theta)^2 (2 * cosd(theta) * sind(theta));
    sind(theta)^2 cosd(theta)^2 (-2 * cosd(theta) * sind(theta));
    (-cosd(theta) * sind(theta)) (cosd(theta) * sind(theta)) ...
    (cosd(theta)^2 - sind(theta)^2)];

% Transformed Reduced compliance matrix ([Sbar])
Sbar = transpose(reducedStressTransformation(theta)) * ...
    S * reducedStressTransformation(theta);

% Reduced strain transformation
reducedStrainTransformation = @(theta) ...
    [cosd(theta)^2 sind(theta)^2 (cosd(theta) * sind(theta));
    sind(theta)^2 cosd(theta)^2 (-cosd(theta) * sind(theta));
    (-2 * cosd(theta) * sind(theta)) (2 * cosd(theta) * sind(theta)) ...
    (cosd(theta)^2 - sind(theta)^2)];

% Transformed Reduced compliance matrix ([Qbar])
Qbar = transpose(reducedStrainTransformation(theta)) * ...
    Q * reducedStrainTransformation(theta);

%% Stress and Strain
if condition == 1
    
    globStress = globSys.stress;
    locStress = reducedStressTransformation(theta) * globStress;
    
    globStrain = Sbar * globStress;
    locStrain = reducedStrainTransformation(theta) * globStrain;
    
elseif condition == 2
    
    locStress = locSys.stress;
    globStress = reducedStressTransformation(-theta) * locSys.stress;
    
    globStrain = Sbar * globStress;
    locStrain = reducedStrainTransformation(theta) * globStrain;
    
elseif condition == 3
    
    globStrain = globSys.strain;
    locStrain = reducedStrainTransformation(theta) * globStrain;
    
    globStress = Qbar * globStrain;
    locStress = reducedStressTransformation(theta) * globStress;
    
elseif condition == 4
    
    locStrain = locSys.strain;
    globStrain = reducedStrainTransformation(-theta) * locStrain;
    
    globStress = Qbar * globStrain;
    locStress = reducedStressTransformation(theta) * globStress;
end


end
