function [laminaStressStrain] = ...
    transformation(globLaminaStress, Sbar, Te, Tsigma)

%% Global and Local Stress and Strain

for ii = 1:size(globLaminaStress, 3)
    
    % Global state of stress
    laminaStressStrain.globStress(:, :, ii) = globLaminaStress(:, :, ii)';
    
    % Local state of stress
    laminaStressStrain.locStress(:, :, ii) = Tsigma(:, :, ii) * ...
        laminaStressStrain.globStress(:, :, ii);
    
    % Global state of strain
    laminaStressStrain.globStrain(:, :, ii) = Sbar(:, :, ii) * ...
        laminaStressStrain.globStress(:, :, ii);
    
    % Local state of strain
    laminaStressStrain.locStrain(:, :, ii) = Te(:, :, ii) * ...
        laminaStressStrain.globStrain(:, :, ii);
    
end
