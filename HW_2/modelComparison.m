function [properties] = modelComparison(Ef, vf, cf, Em, vm)

%% Main Function Body

% Calculating matrix volume ratio
cm = 1 - cf

% Defining all models compared in this subroutine
models = {'Voigt', 'Reuss', 'Hybrid', 'Square Fiber', 'Halpin-Tsai'};

% Pupulating cell array with all model functions
functions = {@voigt, @reuss, @hybrid, @SFM, @halpinTsai};

% Populating the structure with all of the models
for ii = 1:length(models)
    % Assigning model names to structure
    properties(ii).model = models{ii};
    
    % Calling model functions as storing properties in structure
    [properties(ii).E1, ...
     properties(ii).E2, ... 
     properties(ii).v12, ...
     properties(ii).G12] = functions{ii}(Ef, vf, cf, Em, vm , cm);
    
end

    %% Voigt model
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

    %% Reuss Model
    function [E1, E2, v12, G12] = reuss(Ef, vf, cf, Em, vm , cm)
        
        % Function for estimated shear modulus
        G = @(E, v) E / (2 * (1 + v));
        
        % Estimated fiber shear moudulus
        Gf = G(Ef, vf);
        
        % Estimated matrix shear modulus
        Gm = G(Em, vm);
        
        % Compliance matrix for matrix
        S.m = [1 / Em, -vm / Em, 0; ...
               -vm / Em, 1 / Em, 0; ...
               0, 0, 1 / Gm];
        
        % Compliance matrix for fiber
        S.f = [1 / Ef, -vf / Ef, 0; ...
               -vf / Ef, 1 / Ef, 0; ...
               0, 0, 1 / Gf];
           
        % Compliance matrix according to Reuss
        S.r = (cf .* S.f) + (cm .* S.m); 
        
        % Reuss longitudinal modulus
        E1 = 1 / S.r(1, 1);
        
        % Reuss transverse modulus
        E2 = 1 / S.r(2, 2);
        
        % Reuss major Poisson's ratio
        v12 = -S.r(1, 2) * E1;
        
        % Reuss shear modulus
        G12 = 1 / S.r(3, 3);
        
    end

    %% Hybrid model
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

    %% Square Fiber Model
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
    
    %% Halpin-Tsai Model
    function [E1, E2, v12, G12] = halpinTsai(Ef, vf, cf, Em, vm , cm)
        
        % Function for estimated shear modulus
        G = @(E, v) E / (2 * (1 + v));
        
        % Estimated matrix shear modulus
        Gm = G(Em, vm);
        
        % Estimated fiber shear moudulus
        Gf = G(Ef, vf);
        
        % Eta function
        eta = @(Pf, Pm, ep) (Pf - Pm) / (Pf - (ep * Pm));
        
        % Halpin-Tsai E1 = Voigt E1
        E1 = (cf * Ef) + (cm * Em);
        
        % Estimated Parameter for E2-Circular fibers-Square arra7
        estimatedParameter.E2 = 2;
        
        % Eta for E2
        etaE2 = eta(Ef, Em, estimatedParameter.E2);
        
        % Halpin-Tsai E2
        E2 = (Em * (1 + (estimatedParameter.E2 * etaE2 * cf))) / ...
            (1 - (etaE2 * cf));
      
        % Halpin-Tsai v12 = voight v12
        v12 = (cf * vf) + (cm * vm);
        
        % Estimated Parameter for E2-Circular fibers-Square arra
        estimatedParameter.G12 = 1;
        
        % Eta for G12
        etaG12 = eta(Gf, Gm, estimatedParameter.G12);
        
        % Halpin-Tsai G12
        G12 = (Gm * (1 + (estimatedParameter.G12 * etaG12 * cf))) / ...
            (1 - (etaG12 * cf));
         
    end

end
