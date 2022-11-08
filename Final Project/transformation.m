function [locStress, globStress, locStrain, globStrain] = ...
    transformation(globLaminaStress, Qbar, Sbar, Te, Tsigma)

%% Global and Local Stress and Strain

for ii = 1:size(globLaminaStress, 3)
    
    globStress(:, :, ii) = globLaminaStress(:, :, ii)';
    locStress(:, :, ii) = Tsigma(:, :, ii) * globStress(:, :, ii);
    
    globStrain(:, :, ii) = Sbar(:, :, ii) * globStress(:, :, ii);
    locStrain(:, :, ii) = Te(:, :, ii) * globStrain(:, :, ii);
    
end
