function [hygrothermal, laminaStressStrain, superimposedParam, faliure, ...
    effectiveLaminatePropsMech] = main()

%% INPUTS

% % Lamina Properties
% moduli.E1 = 20 * 10^6;
% moduli.E2 = 1.4 * 10^6;
% moduli.G12 = .8 * 10^6;
% moduli.v12 = .3;
% 
% % Strength Properties
% strength.X = 310 * 10 ^ 3; % psi
% strength.Xprime = -310 * 10 ^ 3; % psi
% strength.Y = 9 * 10 ^ 3; % psi
% strength.Yprime = -30 * 10 ^ 3; % psi
% strength.S = 15 * 10 ^ 3; % psi
% strength.Xe = .01555;
% strength.XePrime = -.01555;
% strength.Ye = .006;
% strength.YePrime = -.02;
% strength.Se = .015;

% Lamina Properties
moduli.E1 = 7.4 * 10^6;
moduli.E2 = 1.4 * 10^6;
moduli.G12 = .73 * 10^6;
moduli.v12 = .28;

% Strength Properties
strength.X = 260 * 10 ^ 3; % psi
strength.Xprime = -162 * 10 ^ 3; % psi
strength.Y = 6.2 * 10 ^ 3; % psi
strength.Yprime = -32 * 10 ^ 3; % psi
strength.S = 13.1 * 10 ^ 3; % psi
strength.Xe = .03514;
strength.XePrime = -.02189;
strength.Ye = .00443;
strength.YePrime = -.0229;
strength.Se = .018;

% % Global Hygrothermal Properties
% hygrothermal.alpha = [-.5 * 10 ^ -6, 15 * 10 ^ -6, 0]; % 1/F  

% Global Hygrothermal Properties
hygrothermal.alpha = [3.9 * 10 ^ -6, 16.7 * 10 ^ -6, 0]; % 1/F  

% Temperature conditions
hygrothermal.T0 = 100;
hygrothermal.Tf = 100;

% Loading
loading = [0 0 896 0 0 0]; % 1:3 = forces 4:6 = moments
% Have check in place to ensure this is row vec or col vec

t = .009; % Laminae thickness
ss = [45 -45 0 0 0 0 -45 45]; % Stackup sequence
ssMod = [45 -45 0 0 0 0 0 -45 45];

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

% Checking faliure criteria
[faliure] = faliureCriteria(strength, superimposedParam);

% Effective mechanical and thermal laminate properties
[effectiveLaminatePropsMech] = effectiveLaminateProps(ABD, t, ss, ...
    hygrothermal);

end