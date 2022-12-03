function [globLaminaStress, zMod] = globalLaminaStress(...
    deformationAtMidplane, Qbar, z, ss)

% Initializing matrix for global lamina stress
globLaminaStress = zeros(1, 3, length(ss));

% Modified z-coord mat without 0 position
zMod = z(z ~= 0);

for ii = 1:length(ss)
    
    % Global stress in lamina
    globLaminaStress(:, :, ii) = ...
        Qbar(:, :, ii) * deformationAtMidplane(1:3) + ...
        Qbar(:, :, ii) * zMod(ii) * deformationAtMidplane(4:6);
    
end
