function [properties] = modelComparison(Ef, vf, cf, Em, vm , cm)

% Defining all models compared in this subroutine
models = {'Voigt', 'Reuss', 'Hybrid', 'Square Fiber', 'Halpin-Tsai'};

functions = {@voigt, @reuss, @hybrid, @SFM, @halpinTsai};

% Populating the structure with all of the models
for ii = 1:length(models)
    properties(ii).model = models{ii};
    
    [properties(ii).E1, ...
     properties(ii).E2, ...
     properties(ii).v12, ...
     properties(ii).G12] = functions{ii}(Ef, vf, cf, Em, vm , cm);
    
end

    % Voigt model
    function [E1, E2, v12, G12] = voigt(Ef, vf, cf, Em, vm , cm)
        
        % Function for estimated shear modulus
        G = @(E, v) E / (2 * (1 + v));
        
        % Estimated fiber shear moudulus
        Gf = G(Ef, vf);
        
        % Estimated matrix shear modulus
        Gm = G(Em, vm);
        
        % Stifness matrix for matrix
        C.m = [Em / (1 - (vm ^ 2)), (vm * Em) / (1 - (vm ^ 2)), ...
              0; (vm * Em) / (1 - (vm ^ 2)), Em / (1 - (vm ^ 2)), ...
              0; 0, 0, Gm];
          
        % Stiffness matrix for fibers
        C.f = [Ef / (1 - (vf ^ 2)), (vf * Ef) / (1 - (vf ^ 2)), ...
              0; (vf * Ef) / (1 - (vf ^ 2)), Ef / (1 - (vf ^ 2)), ...
              0; 0, 0, Gf];
        
        % Stiffness matrix according to Voigt
        C.v = (cf .* C.f) + (cm .* C.m);
        
        % Compliance matrix according to Voigt
        S.v = inv(C.v);
        
        % Voigt longitudinal modulus
        E1 = 1 / S.v(1, 1);
        
        % Voigt transverse modulus
        E2 = 1 / S.v(2, 2);
        
        % Voigt major Poisson's ratio
        v12 = -S.v(2, 1) * E1;
        
        % Voigt shear modulus
        G12 = 1 /S.v(3, 3);
        
    end

    function [E1, E2, v12, G12] = reuss(Ef, vf, cf, Em, vm , cm)
        
        E1 = 1;
        E2 = 1;
        v12 = 1;
        G12 = 1;
        
    end

    % Hybrid model
    function [E1, E2, v12, G12] = hybrid(Ef, vf, cf, Em, vm , cm)
        
        % Hybirid longitudinal modulus
        E1 = (cf * Ef) + (cm * Em);
        
        % Hybrid transverse modulus
        E2 = 1 / (((cf / Ef) + (cm / Em)) - (((cf * cm) / (Ef * Em)) * ...
            ((((vf * Em) - (vm * Em)) ^ 2) / ((cm * Em) + (cf * Ef)))));
        
        % Function for estimated shear modulus
        G = @(E, v) E / (2 * (1 + v));
        
        % Estimated fiber shear moudulus
        Gf = G(Ef, vf);
        
        % Estimated matrix shear modulus
        Gm = G(Em, vm);
        
        % Hybrid in-palne shear modulus
        G12 = 1 / ((cf * (1 / Gf)) + (cm * (1 / Gm)));
        
        % Hybrid major Poisson's ratio
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
        
        E1 = 1;
        E2 = 1;
        v12 = 1;
        G12 = 1;
        
    end

end
