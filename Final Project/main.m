function [hygrothermal, laminaStressStrain, superimposedParam] = main()

%% INPUTS

% Constituent Properties
E1f = 73;
E2f = 20;
vf = .22;
cf = .55;

Em = 3.5;
vm = .35;

% Global Hygrothermal Properties
hygrothermal.alpha = [-.5 * 10 ^ -6, 15 * 10 ^ -6, 0]; % 1/F  

% Temperature conditions
hygrothermal.T0 = 350;
hygrothermal.Tf = 100;

% Moisture conditions
hygrothermal.C0 = 0;
hygrothermal.Cf = 0;

% Loading
loading = [100 100 100 200 100 300]; % 1:3 = forces 4:6 = moments
% Have check in place to ensure this is row vec or col vec

t = .1; % Laminae thickness
ss = [0 30 -30 -30 30 30]; % Stackup sequence

%% Main Code Body

% Calculating lamina properties from micromechanical properties
% Utilizing SFM
[moduli] = SFM(E1f, E2f, vf, cf, Em, vm);

% Create Lamina strct with Qbar
[Qbar, S, Tepsilon] = reducedTransformedStiffnessMat(moduli, ss);

% Add Sbar to lamina struct
[Sbar, Tsigma] = reducedTransformedComplianceMat(S, ss);

% Create laminate strcture with deformationAtMidplane & ABD matrix
[deformationAtMidplane, z, ABD] = midplaneDeformation(loading, Qbar, ss, t);

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



end