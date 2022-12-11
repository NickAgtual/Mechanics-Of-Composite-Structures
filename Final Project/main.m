function [hygrothermal, laminaStressStrain, superimposedParam, Qbar] = main()

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
ssMod = [0 30 -30 -30 -30 30 0];

%% Main Code Body

% Create Lamina strct with Qbar
[Qbar, S, Tepsilon] = reducedTransformedStiffnessMat(moduli, ssMod);

% Add Sbar to lamina struct
[Sbar, Tsigma] = reducedTransformedComplianceMat(S, ssMod);

% Create laminate strcture with deformationAtMidplane & ABD matrix
[deformationAtMidplane, z, ABD] = midplaneDeformation(loading, Qbar, ...
    ssMod, t, ss);

% Add global stress to lamina struct
[globLaminaStress] = globalLaminaStress(deformationAtMidplane, ...
    Qbar, z, ssMod);

% Add all remaining stress & strain to lamina struct
[laminaStressStrain] = transformation(globLaminaStress, Sbar, Tepsilon, ...
    Tsigma);

% Calculating hygrothermal stresses
[hygrothermal] = hygrothermalEffetcs(ssMod, hygrothermal, z, Tepsilon, ...
    Tsigma, ABD, Qbar, ss);

% Superimposing stress and strain
[superimposedParam] = superposition(laminaStressStrain, hygrothermal);

% Plotting global stress and strain
stressStrainPlots(superimposedParam, z)




end