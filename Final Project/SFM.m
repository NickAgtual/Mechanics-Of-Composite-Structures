function [moduli] = SFM(E1f, E2f, vf, cf, Em, vm)

% Modified volume ratios
cA = sqrt(cf);
cB = 1 - cA;

% Function for estimated shear modulus
G = @(E, v) E / (2 * (1 + v));

% Estimated matrix shear modulus
Gm = G(Em, vm);

% Estimated fiber shear moudulus
G12f = G(E1f, vf);

% Plain reduced stiffness matrix
Q.m = [Em / (1 - (vm ^ 2)), (vm * Em) / (1 - (vm ^ 2)), ...
    0; (vm * Em) / (1 - (vm ^ 2)), Em / (1 - (vm ^ 2)), ...
    0; 0, 0, Gm];

% Hybrid in-plane modulus
hybrid.E1 = (cA * E1f) + (cB * Em);

% Hybrid major Poisson's ratio
hybrid.v12 = (cA * vf) + (cB * vm);

% Hybrid shear modulus
hybrid.G12 = 1 / ((cA / G12f) + (cB / Gm));

% Hybrid transverse modulus
hybrid.E2 = (cA / E2f) + ((cB / Em) * (1 - (vm ^ 2)));

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
moduli.E1 = ((Q.SFM(1, 1) * Q.SFM(2, 2)) - (Q.SFM(1, 2) ^ 2)) / ...
    Q.SFM(2, 2);

% SFM transverse modulus
moduli.E2 = ((Q.SFM(1, 1) * Q.SFM(2, 2)) - (Q.SFM(1, 2) ^ 2)) / ...
    Q.SFM(1, 1);

% SFM shear modulus
moduli.G12 = Q.SFM(3, 3);

% SFM major Poisson's ratio
moduli.v12 = Q.SFM(1, 2) / Q.SFM(2, 2);



end