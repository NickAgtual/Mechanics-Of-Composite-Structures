function [globLaminaStress] = globalLaminaStress(...
    deformationAtMidplane, Qbar, t, z)

globLaminaStress = zeros(1, 3, length(t));

for ii = 1:length(t)
    
    globLaminaStress(:, :, ii) = ...
        Qbar(:, :, ii) * deformationAtMidplane(1:3) + ...
        Qbar(:, :, ii) * z(ii) * deformationAtMidplane(4:6);

end
