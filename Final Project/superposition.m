function [superimposedParam] = superposition(laminaStressStrain, ...
    hygrothermal)

% structure fields to loop over
structFields = {'globStress', 'locStress', 'globStrain', 'locStrain'};

% Superimposing loc/global stress/strain due to mechanical and hygrothermal
% loading
for ii = 1:length(structFields)
    
    superimposedParam.(structFields{ii}) = ...
        hygrothermal.(structFields{ii}) + ...
        laminaStressStrain.(structFields{ii});
    
end

end

