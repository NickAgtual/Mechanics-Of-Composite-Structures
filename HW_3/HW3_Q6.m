clear; clc; close all

%% Material Properties

% Longitudinal modulus
lamina.E1 = 147 * 10 ^ 9; % Pa

% Transverse moudulus
lamina.E2 = 10.3 * 10 ^ 9; % Pa

% In-plane shear modulus
lamina.G12 = 7 * 10 ^ 9; % Pa

% Poisson's ratio
lamina.v12 = .27;

% Stacking sequences
laminate.layup = [30 -30 -30 30];

% Thickness of lamina
lamina.thickness = .127; % mm

hygrothermal.alpha1 = -.9 * 10 ^ -6;
hygrothermal.alpha2 = 27 * 10 ^ -6;
hygrothermal.alpha = [hygrothermal.alpha1 hygrothermal.alpha2 0];

% Temperature difference
hygrothermal.deltaT = 56;

% Lamina coordinates within laminate

z = zeros(1, length(laminate.layup) / 2);

for ii = 1: length(laminate.layup) / 2
    
    z(ii) = lamina.thickness * ii;
    
end

laminate.zCoord = [-flip(z) 0 z];
laminate.zCoordMod = laminate.zCoord(laminate.zCoord ~= 0);

%% 1) Calculate the A, B, D, a, b, d Matrices for the Layup [45/90/30/0]

% Compliance matrix
matrices.S = [(1 / lamina.E1) (-lamina.v12 / lamina.E1) 0;
              (-lamina.v12 / lamina.E1) (1 / lamina.E2) 0;
              0 0 (1 / lamina.G12)];
          
% Stiffness matrix
matrices.Q = inv(matrices.S);

% Initializing Qbar matrix
matrices.Qbar = zeros(3, 3, length(laminate.layup));

% Initializing strain transformation matrix
matrices.Tepsilon = zeros(3, 3, length(laminate.layup));

for ii = 1:length(laminate.layup)
    
    % Strain transformation matrix
    matrices.Tepsilon(:, :, ii) = [cosd(laminate.layup(ii))^2 ...
        sind(laminate.layup(ii))^2 ...
        cosd(laminate.layup(ii)) * sind(laminate.layup(ii));
        sind(laminate.layup(ii))^2 ...
        cosd(laminate.layup(ii))^2 ...
        -cosd(laminate.layup(ii)) * sind(laminate.layup(ii));
        -2 * cosd(laminate.layup(ii)) * sind(laminate.layup(ii)) ...
        2 * cosd(laminate.layup(ii)) * sind(laminate.layup(ii)) ...
        cosd(laminate.layup(ii))^2 - sind(laminate.layup(ii))^2];
    
    % Reduced transformed stiffness matrix
    matrices.Qbar(:, :, ii) = transpose(matrices.Tepsilon(:, :, ii)) * ...
        matrices.Q * matrices.Tepsilon(:, :, ii);
                    
end

%% Global Hygrothermal Properties & Free Strain
for ii = 1:length(laminate.layup)
    
    % Global hygrothermal properties
    hygrothermal.alphaGlobal(:, :, ii) = inv(matrices.Tepsilon(:, :, ii)) * ...
        hygrothermal.alpha';
    
    % Free strain
    hygrothermal.freeStrain(:, :, ii) = ...
        (hygrothermal.alphaGlobal(:, :, ii) .* hygrothermal.deltaT);
    
end

%% Hygrothermal Loads

% Initializing hygrothermal loads
[hygrothermal.Nprelim, hygrothermal.Mprelim] = deal(zeros(3, 1));

% Calculating forces and moments due to hygrothermal conditions
for ii = 1:length(laminate.layup)

    % Forces (does not include mutliplication of temp diff.)
    hygrothermal.Nprelim = hygrothermal.Nprelim + ((matrices.Qbar(:, :, ii) ...
        * hygrothermal.alphaGlobal(:, :, ii) .* ...
        (laminate.zCoord(ii+1) - laminate.zCoord(ii)))); 
    
    % Moments (does not include mutliplication of temp diff.)
    hygrothermal.Mprelim = hygrothermal.Mprelim + ((matrices.Qbar(:, :, ii) * ...
        hygrothermal.alphaGlobal(:, :, ii) * ((laminate.zCoord(ii + 1) ^ 2) ...
        - (laminate.zCoord(ii) ^ 2))));
end

% Thermal Loads
hygrothermal.N = hygrothermal.Nprelim * hygrothermal.deltaT;
hygrothermal.M = hygrothermal.Mprelim * hygrothermal.deltaT * .5;

for ii = 1:length(laminate.layup)
    
    % Strain transformation matrix
    matrices.Tepsilon(:, :, ii) = [cosd(laminate.layup(ii))^2 ...
        sind(laminate.layup(ii))^2 ...
        cosd(laminate.layup(ii)) * sind(laminate.layup(ii));
        sind(laminate.layup(ii))^2 ...
        cosd(laminate.layup(ii))^2 ...
        -cosd(laminate.layup(ii)) * sind(laminate.layup(ii));
        -2 * cosd(laminate.layup(ii)) * sind(laminate.layup(ii)) ...
        2 * cosd(laminate.layup(ii)) * sind(laminate.layup(ii)) ...
        cosd(laminate.layup(ii))^2 - sind(laminate.layup(ii))^2];
    
    % Reduced transformed stiffness matrix
    matrices.Qbar(:, :, ii) = transpose(matrices.Tepsilon(:, :, ii)) * ...
        matrices.Q * matrices.Tepsilon(:, :, ii);
                    
end

% Initializing ABD matrices
[matrices.A, matrices.B, matrices.D] = deal(zeros(3, 3));

% Solving for ABD matrices
for ii = 1:3
    for jj = 1:3
        for kk = 2:length(laminate.layup) + 1
            
            matrices.A(ii, jj) = matrices.A(ii, jj) + ...
                (matrices.Qbar(ii, jj, kk-1) * ...
                (laminate.zCoord(kk) - laminate.zCoord(kk-1)));
            
            matrices.B(ii, jj) = matrices.B(ii, jj) + ...
                (.5 * (matrices.Qbar(ii, jj, kk-1) ...
                * (laminate.zCoord(kk)^2 - laminate.zCoord(kk-1)^2)));
            
            matrices.D(ii, jj) = matrices.D(ii, jj) + ((1/3) * ...
                (matrices.Qbar(ii, jj, kk-1)...
                 * (laminate.zCoord(kk)^3 - laminate.zCoord(kk-1)^3)));
            
        end
    end
end

% ABD matrix
matrices.ABD = [matrices.A matrices.B; matrices.B matrices.D];

% abd matrix
matrices.abd = inv(matrices.ABD);

% Hygrothermal midplain deformation
% Since this only accounts for temperature change, the total midplane def. 
% is equivalent to the midplain strain and curverature due to hygrothermal
% conditions
hygrothermal.midplaneDeformation = inv(matrices.ABD) * ...
    [hygrothermal.N; hygrothermal.M];

for ii = 1:length(laminate.layup)
    
    hygrothermal.globStrain(:, :, ii) = ...
        hygrothermal.midplaneDeformation(1:3) + ...
        (laminate.zCoordMod(ii) * hygrothermal.midplaneDeformation(4:end)) - ...
        hygrothermal.freeStrain(:, :, ii);
end

% Hygrothermal stress (lamina) and local strain
for ii = 1:length(laminate.layup)

    % Global Stress
    hygrothermal.globStress(:, :, ii) = matrices.Qbar(:, :, ii) * ...
        (hygrothermal.midplaneDeformation(1:3) + (laminate.zCoordMod(ii) * ...
        hygrothermal.midplaneDeformation(4:end)) - ...
        hygrothermal.freeStrain(:, :, ii));
    
    % Local Stress
    hygrothermal.locStress(:, :, ii) = Tsigma(:, :, ii) * ...
        hygrothermal.globStress(:, :, ii);
    
    % Local strain
    hygrothermal.locStrain(:, :, ii) = Tepsilon(:, :, ii) * ...
        hygrothermal.globStrain(:, :, ii);
    
end
