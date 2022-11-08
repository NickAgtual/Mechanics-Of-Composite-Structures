function [deformationAtMidplane, globLaminaStress] = main()

%% INPUTS

% Constituent Properties
E1f = 73;
E2f = 20;
vf = .22;
cf = .55;

Em = 3.5;
vm = .35;

% Loading
loading = [100 100 100 200 100 300]; % 1:3 = forces 4:6 = moments
% Have check in place to ensure this is row vec or col vec
t = [.1 .1 .1]; % Laminae thickness
ss = [0 0 90 90]; % Stackup sequence

%% Main Code Body

[moduli] = SFM(E1f, E2f, vf, cf, Em, vm);

[Qbar, S] = reducedTransformedStiffnessMat(moduli, ss);

[Sbar] = reducedTransformedComplianceMat(S, ss);

[deformationAtMidplane, z] = midplaneDeformation(loading, Qbar, t);

[globLaminaStress] = globalLaminaStress(deformationAtMidplane, Qbar, t, z);



end