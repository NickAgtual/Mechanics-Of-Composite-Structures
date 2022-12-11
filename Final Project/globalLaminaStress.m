function [globLaminaStress] = globalLaminaStress(...
    deformationAtMidplane, Qbar, z, ssMod)

% Initializing matrix for global lamina stress
globLaminaStress = zeros(1, 3, length(ssMod));

for ii = 1:length(ssMod)
    
    % Global stress in lamina
    globLaminaStress(:, :, ii) = ...
        Qbar(:, :, ii) * deformationAtMidplane(1:3) + ...
        Qbar(:, :, ii) * z(ii) * deformationAtMidplane(4:6);
    
end
