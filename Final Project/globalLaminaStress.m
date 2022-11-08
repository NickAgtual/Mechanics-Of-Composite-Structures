function [globLaminaStress] = globalLaminaStress(...
    deformationAtMidplane, Qbar, t, z)

globLaminaStress = zeros(1, 3, ii);

for ii = 1:length(t)
    
    globLaminaStress(:, :, ii) = ...
        Qbar(:, :, ii) * deformationAtMidplane(1:3) + ...
        Qbar(:. :, ii) * z * deformationAtMidplane(4:6);

end
