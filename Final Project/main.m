function [hygrothermal] = main()

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
hygrothermal.beta = [0 0 0];  

% Temperature conditions
hygrothermal.T0 = 350;
hygrothermal.Tf = 100;

% Moisture conditions
hygrothermal.C0 = 0;
hygrothermal.Cf = 0;

% Loading
loading = [100 100 100 200 100 300]; % 1:3 = forces 4:6 = moments
% Have check in place to ensure this is row vec or col vec

t = [.1 .1 .1]; % Laminae thickness
ss = [0 0 90 90]; % Stackup sequence

%% Main Code Body

[moduli] = SFM(E1f, E2f, vf, cf, Em, vm);

% Create Lamina strct with Qbar
[Qbar, S, Te] = reducedTransformedStiffnessMat(moduli, ss);

% Add Sbar to lamina struct
[Sbar, Tsigma] = reducedTransformedComplianceMat(S, ss);

% Create laminate strcture with deformationAtMidplane & ABD matrix
[deformationAtMidplane, z, ABD] = midplaneDeformation(loading, Qbar, t);

% Calculating hygrothermal stresses
[hygrothermal] = hygrothermalEffetcs(ss, hygrothermal, z, Te, ABD, Qbar)

% Add global stress to lamina struct
[globLaminaStress] = globalLaminaStress(deformationAtMidplane, Qbar, t, z);

% Add all remaining stress & strain to lamina struct
[laminaStressStrain] = transformation(globLaminaStress, Sbar, Te, Tsigma);



end