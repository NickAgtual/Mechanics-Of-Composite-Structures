function [properties] = modelComparison(Ef, vf, cf, Em, vm , cm)

% Defining all models compared in this subroutine
models = {'Voigt', 'Reuss', 'Hybrid', 'Square Fiber', 'Halpin-Tsai'};

% Populating the structure with all of the models
for ii = 1:length(models)
    properties(ii).model = models{ii};
    
end


    % Voigt model
    function [E1, E2, G12, v12] = voigt(Ef, vf, cf, Em, vm , cm)
        
        x = 1;
     
    end

    function [E1, E2, G12, v12] = reuss(Ef, vf, cf, Em, vm , cm)
        
        x = 1;
        
    end

    % Hybrid model
    function [E1, E2, G12, v12] = hybrid(Ef, vf, cf, Em, vm , cm)
        
        % Longitudinal modulus
        E1 = (cf * Ef) * (cm * Em);
        
        % Transverse modulus
        E2 = (cf / Ef) + ((cm / Em) * (1 - (vm ^ 2)));
        
        % Function for estimated shear modulus
        G = @(E, v) E / (2 * (1 + v));
        
        % Estimated fiber shear moudulus
        Gf = G(Ef, vf);
        
        % Estimated matrix shear modulus
        Gm = G(Em, vm);
        
        % In-palne shear modulus
        G12 = (cf * (1 / Gf)) + (cm * (1 / Gm));
        
        % Major Poisson's ratio
        v12 = (cf * vf) + (cm * vm);
        
    end

    function [E1, E2, G12, v12] = sfm(Ef, vf, cf, Em, vm , cm)
        
        % Function for estimated shear modulus
        G = @(E, v) E / (2 * (1 + v));
        
        % Estimated matrix shear modulus
        Gm = G(Em, vm);
        
        % Plain reduced stiffness matrix
        Q.m = [Em / (1 - (vm ^ 2)), (vm * Em) / (1 - (vm ^ 2)), ...
              0; (vm * Em) / (1 - (vm ^ 2)), Em / (1 - (vm ^ 2)), ...
              0; 0, 0, Gm];
        
        % Hybrid reduced stiffness matrix
        Q.h = 1;
        
    end

    function [E1, E2, G12, v12] = halpinTsai(Ef, vf, cf, Em, vm , cm)
        
        x = 1;
        
    end

end
