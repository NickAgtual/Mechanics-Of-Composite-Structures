function [properties] = modelComparison(Ef, vf, cf, Em, vm , cm)

% Defining all models compared in this subroutine
models = {'Voigt', 'Reuss', 'Hybrid', 'Square Fiber', 'Halpin-Tsai'};

% Populating the structure with all of the models
for ii = 1:length(models)
    properties(ii).model = models{ii};
    
end

for ii = 1:length(models)
    switch ii
        case 1
            voigt(Ef, vf, cf, Em, vm , cm)
        case 2
            z = 1;
        case 3
            [properties(ii).E1, ...
             properties(ii).E2, ...
             properties(ii).v12, ...
             properties(ii).G12] = hybrid(Ef, vf, cf, Em, vm, cm);
        case 4
            [properties(ii).E1, ...
             properties(ii).E2, ...
             properties(ii).v12, ...
             properties(ii).G12] = SFM(Ef, vf, cf, Em, vm, cm);
        case 5
            z = 1;
    end
end

    % Voigt model
    function [E1, E2, v12, G12] = voigt(Ef, vf, cf, Em, vm , cm)
        
                % Plain reduced stiffness matrix
        C = [Em / (1 - (vm ^ 2)), (vm * Em) / (1 - (vm ^ 2)), ...
              0; (vm * Em) / (1 - (vm ^ 2)), Em / (1 - (vm ^ 2)), ...
              0; 0, 0, Gm];
          
        S = inv(C)
        
        
     
    end

    function [E1, E2, v12, G12] = reuss(Ef, vf, cf, Em, vm , cm)
        
        x = 1;
        
    end

    % Hybrid model
    function [E1, E2, v12, G12] = hybrid(Ef, vf, cf, Em, vm , cm)
        
        % Longitudinal modulus
        E1 = (cf * Ef) + (cm * Em);
        
        % Transverse modulus
        E2 = 1 / (((cf / Ef) + (cm / Em)) - (((cf * cm) / (Ef * Em)) * ...
            ((((vf * Em) - (vm * Em)) ^ 2) / ((cm * Em) + (cf * Ef)))));
        
        % Function for estimated shear modulus
        G = @(E, v) E / (2 * (1 + v));
        
        % Estimated fiber shear moudulus
        Gf = G(Ef, vf);
        
        % Estimated matrix shear modulus
        Gm = G(Em, vm);
        
        % In-palne shear modulus
        G12 = 1 / ((cf * (1 / Gf)) + (cm * (1 / Gm)));
        
        % Major Poisson's ratio
        v12 = (cf * vf) + (cm * vm);
        
    end

    function [E1, E2, v12, G12] = SFM(Ef, vf, cf, Em, vm , cm)
        
        % Modified volume ratios
        cA = sqrt(cf);
        cB = 1 - cA;
        
        % Function for estimated shear modulus
        G = @(E, v) E / (2 * (1 + v));
        
        % Estimated matrix shear modulus
        Gm = G(Em, vm);
        
        % Estimated fiber shear moudulus
        Gf = G(Ef, vf);
        
        % Plain reduced stiffness matrix
        Q.m = [Em / (1 - (vm ^ 2)), (vm * Em) / (1 - (vm ^ 2)), ...
              0; (vm * Em) / (1 - (vm ^ 2)), Em / (1 - (vm ^ 2)), ...
              0; 0, 0, Gm];
          
        % Hybrid in-plane modulus
        hybrid.E1 = (cA * Ef) + (cB * Em);
        
        % Hybrid major Poisson's ratio
        hybrid.v12 = (cA * vf) + (cB * vm);
        
        % Hybrid shear modulus
        hybrid.G12 = 1 / ((cA / Gf) + (cB / Gm));
        
        % Hybrid transverse modulus
        hybrid.E2 = (cA / Ef) + ((cB / Em) * (1 - (vm ^ 2)));
        
        % Hybrid minor Poisson's ratio 
        hybrid.v21 = (hybrid.v12 * hybrid.E2) / hybrid.E1;
        
        % Hybrid reduced stiffness matrix
        Q.h = [hybrid.E1 / (1 - (hybrid.v12 * hybrid.v21)), ...
            (hybrid.v12 * hybrid.E2) / (1 - (hybrid.v12 * hybrid.v21)), ...
            0; (hybrid.v12 * hybrid.E2) / ...
            (1 - (hybrid.v12 * hybrid.v21)), hybrid.E2 / ...
            (1 - (hybrid.v12 * hybrid.v21)), 0; 0, 0, hybrid.G12];
        
        % Reduced stiffness matrix according to SFM
        Q.SFM = (cA .* Q.h) + (cB .* Q.m);
        
        % SFM longitudinal modulus
        E1 = ((Q.SFM(1, 1) * Q.SFM(2, 2)) - (Q.SFM(1, 2) ^ 2)) / ...
            Q.SFM(2, 2);
        
        % SFM transverse modulus
        E2 = ((Q.SFM(1, 1) * Q.SFM(2, 2)) - (Q.SFM(1, 2) ^ 2)) / ...
            Q.SFM(1, 1);
        
        % SFM shear modulus
        G12 = Q.SFM(3, 3);
        
        % SFM major Poisson's ratio
        v12 = Q.SFM(1, 2) / Q.SFM(2, 2);
    end

    function [E1, E2, G12, v12] = halpinTsai(Ef, vf, cf, Em, vm , cm)
        
        x = 1;
        
    end

end
