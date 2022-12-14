clear; close all; clc

%% Contstituent Material Properties
Ef = 395;
vf = .2;
Em = 3.7;
vm = .35;

%% Main Function Body

% Creating range for fiber volume ratio
cf = linspace(0, 1, 100);

% Calculating matrix volume ratio
cm = 1 - cf;

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

%% Plots

% Field names
fn = fieldnames(properties);

for ii = 1:length(fn) - 1
    for jj = 1: length(models)
        
        % Creating new figure
        figure(ii)
        
        % Plotting
        plot(cf, properties(jj).(fn{ii + 1}), 'DisplayName', ...
            properties(jj).model)
        
        % Plot properties
        hold on
        grid on
        grid minor
        
        % Plot Descriptors
        titleText = [fn{ii + 1}, ' vs. Fiber Volume Ratio'];
        title(titleText, 'fontsize', ...
            16, 'Interpreter', 'Latex')
        xlabel('$c_f$', 'fontsize', 14, 'Interpreter', 'Latex')
        ylabel(fn{ii + 1},...
            'fontsize', 14, 'Interpreter', 'Latex')
        legend('location', 'best', 'Interpreter', 'Latex')
        
    end
end

    %% Voigt model
    function [E1, E2, v12, G12] = voigt(Ef, vf, cf, Em, vm , cm)
        
         % Initialzing all modulus vectors
        [E1, E2, G12, v12] = deal(zeros(1, length(cf)));
    
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
          
        for jj = 1:length(cf)
        
            % Stiffness matrix according to Voigt
            C.v = (cf(jj) .* C.f) + (cm(jj) .* C.m);

            % Compliance matrix according to Voigt
            S.v = inv(C.v);

            % Voigt longitudinal modulus
            E1(jj) = 1 / S.v(1, 1);

            % Voigt transverse modulus
            E2(jj) = 1 / S.v(2, 2);

            % Voigt major Poisson's ratio
            v12(jj) = -S.v(2, 1) * E1(jj);

            % Voigt shear modulus
            G12(jj) = 1 /S.v(3, 3);
        
        end
        
    end

    %% Reuss Model
    function [E1, E2, v12, G12] = reuss(Ef, vf, cf, Em, vm , cm)
    
         % Initialzing all modulus vectors
        [E1, E2, G12, v12] = deal(zeros(1, length(cf)));
        
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
        
        for jj = 1:length(cf)
            
            % Compliance matrix according to Reuss
            S.r = (cf(jj) .* S.f) + (cm(jj) .* S.m); 

            % Reuss longitudinal modulus
            E1(jj) = 1 / S.r(1, 1);

            % Reuss transverse modulus
            E2(jj) = 1 / S.r(2, 2);

            % Reuss major Poisson's ratio
            v12(jj) = -S.r(1, 2) * E1(jj);

            % Reuss shear modulus
            G12(jj) = 1 / S.r(3, 3);

        end
        
    end

    %% Hybrid model
    function [E1, E2, v12, G12] = hybrid(Ef, vf, cf, Em, vm , cm)
    
        % Initialzing all modulus vectors
        [E1, E2, G12, v12] = deal(zeros(1, length(cf)));
        
        for jj = 1:length(cf)
        
            % Hybirid longitudinal modulus
            E1(jj) = (cf(jj) * Ef) + (cm(jj) * Em);
            
            % Hybrid transverse modulus
            E2(jj) = 1 / (((cf(jj) / Ef) + (cm(jj) / Em)) - ...
                (((cf(jj) * cm(jj)) / (Ef * Em)) * ...
                ((((vf * Em) - (vm * Em)) ^ 2) / ...
                ((cm(jj) * Em) + (cf(jj) * Ef)))));
            
            % Function for estimated shear modulus
            G = @(E, v) E / (2 * (1 + v));
            
            % Estimated fiber shear moudulus
            Gf = G(Ef, vf);
            
            % Estimated matrix shear modulus
            Gm = G(Em, vm);
            
            % Hybrid in-palne shear modulus
            G12(jj) = 1 / ((cf(jj) * (1 / Gf)) + (cm(jj) * (1 / Gm)));
            
            % Hybrid major Poisson's ratio
            v12(jj) = (cf(jj) * vf) + (cm(jj) * vm);
        
        end
        
    end

    %% Square Fiber Model
    function [E1, E2, v12, G12] = SFM(Ef, vf, cf, Em, vm , ~)
    
    	% Initialzing all modulus vectors
        [E1, E2, G12, v12] = deal(zeros(1, length(cf)));
    
        for jj = 1:length(cf)
        
        % Modified volume ratios
        cA = sqrt(cf(jj));
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
        E1(jj) = ((Q.SFM(1, 1) * Q.SFM(2, 2)) - (Q.SFM(1, 2) ^ 2)) / ...
            Q.SFM(2, 2);
        
        % SFM transverse modulus
        E2(jj) = ((Q.SFM(1, 1) * Q.SFM(2, 2)) - (Q.SFM(1, 2) ^ 2)) / ...
            Q.SFM(1, 1);
        
        % SFM shear modulus
        G12(jj) = Q.SFM(3, 3);
        
        % SFM major Poisson's ratio
        v12(jj) = Q.SFM(1, 2) / Q.SFM(2, 2);
        
        end
        
    end
    
    %% Halpin-Tsai Model
    function [E1, E2, v12, G12] = halpinTsai(Ef, vf, cf, Em, vm , cm)
    
         % Initialzing all modulus vectors
        [E1, E2, G12, v12] = deal(zeros(1, length(cf)));
        
        % Function for estimated shear modulus
        G = @(E, v) E / (2 * (1 + v));
        
        % Estimated matrix shear modulus
        Gm = G(Em, vm);
        
        % Estimated fiber shear moudulus
        Gf = G(Ef, vf);
        
        % Eta function
        eta = @(Pf, Pm, ep) (Pf - Pm) / (Pf + (ep * Pm));
        
        for ii = 1:length(cf)
            
            % Halpin-Tsai E1 = Voigt E1
            E1(ii) = (cf(ii) * Ef) + (cm(ii) * Em);
            
            % Estimated Parameter for E2-Circular fibers-Square arra7
            estimatedParameter.E2 = 2;
            
            % Eta for E2
            etaE2 = eta(Ef, Em, estimatedParameter.E2);
            
            % Halpin-Tsai E2
            E2(ii) = (Em * (1 + (estimatedParameter.E2 * etaE2 * ...
                cf(ii)))) / (1 - (etaE2 * cf(ii)));
            
            % Halpin-Tsai v12 = voight v12
            v12(ii) = (cf(ii) * vf) + (cm(ii) * vm);
            
            % Estimated Parameter for E2-Circular fibers-Square array
            estimatedParameter.G12 = 1;
            
            % Eta for G12
            etaG12 = eta(Gf, Gm, estimatedParameter.G12);
            
            % Halpin-Tsai G12
            G12(ii) = (Gm * (1 + (estimatedParameter.G12 * etaG12 .* ...
                cf(ii)))) /(1 - (etaG12 * cf(ii)));
            
        end
        
    end
   



