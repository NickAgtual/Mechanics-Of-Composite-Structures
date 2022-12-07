function [hygrothermal, laminaStressStrain, superimposedParam] = main()

%% INPUTS

% Constituent Properties
moduli.E1 = 20 * 10^6;
moduli.E2 = 1.4 * 10^6;
moduli.G12 = .8 * 10^6;
moduli.v12 = .3;

% Global Hygrothermal Properties
hygrothermal.alpha = [-.5 * 10 ^ -6, 15 * 10 ^ -6, 0]; % 1/F  

% Temperature conditions
hygrothermal.T0 = 0;
hygrothermal.Tf = -250;

% Loading
loading = [0 10 0 8 -2 3]; % 1:3 = forces 4:6 = moments
% Have check in place to ensure this is row vec or col vec

t = .005; % Laminae thickness
ss = [0 30 -30 -30 30 0]; % Stackup sequence

%% Main Code Body

% Create Lamina strct with Qbar
[Qbar, S, Tepsilon] = reducedTransformedStiffnessMat(moduli, ss);

% Add Sbar to lamina struct
[Sbar, Tsigma] = reducedTransformedComplianceMat(S, ss);

% Create laminate strcture with deformationAtMidplane & ABD matrix
[deformationAtMidplane, z, ABD] = midplaneDeformation(loading, Qbar, ...
    ss, t);

% Add global stress to lamina struct
[globLaminaStress, zMod] = globalLaminaStress(deformationAtMidplane, ...
    Qbar,z, ss);

% Add all remaining stress & strain to lamina struct
[laminaStressStrain] = transformation(globLaminaStress, Sbar, Tepsilon, ...
    Tsigma);

% Calculating hygrothermal stresses
[hygrothermal] = hygrothermalEffetcs(ss, hygrothermal, z, Tepsilon, ...
    Tsigma, ABD, Qbar, zMod);

% Superimposing stress and strain
[superimposedParam] = superposition(laminaStressStrain, hygrothermal);

% Plotting global stress and strain
stressStrainPlots(superimposedParam, zMod)


end