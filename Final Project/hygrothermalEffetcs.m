function [hygrothermal] = hygrothermalEffetcs(ss, hygrothermal, z, Te)

%% Local Hygrothermal Properties
for ii = 1:length(ss)
    
    hygrothermal.alphaLocal(:, :, ii) = Te(:, :, ii) * hygrothermal.alpha';
    hygrothermal.betaLocal(:, :, ii) = Te(:, :, ii) * hygrothermal.beta';
   
end
