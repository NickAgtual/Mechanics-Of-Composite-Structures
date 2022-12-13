%% HW3_Q7
clear; clc; close all
%% Given Properties and Conditions
% Hygrothermal conditions
hygrothermal.T1 = 180; % deg C
hygrothermal.T2 = 30; % deg C
hygrothermal.deltaT = hygrothermal.T2 - hygrothermal.T1;

% Mechanical properties
lamina.E1 = 140; % GPa
lamina.E2 = 10; % GPa
lamina.G12 = 6; % GPa
lamina.v12 = .34;
lamina.t = 1 * 10 ^ -3; % m

% Hygrothermal properties
hygrothermal.alpha1 = 0;
hygrothermal.beta1 = 0;
hygrothermal.alpha2 = 30 * 10 ^ -6; % 1 / deg C
hygrothermal.beta2 = .55;

% Combining hygrothermal properties
hygrothermal.alpha = [hygrothermal.alpha1 hygrothermal.alpha2 0];
hygrothermal.beta = [hygrothermal.beta1 hygrothermal.beta2 0];

layup = [0 90 -90 0];
laminate.zCoord = [-2*lamina.t -lamina.t 0 lamina.t 2*lamina.t];
zMod = laminate.zCoord(laminate.zCoord ~= 0);

%% Reduced Transfromed Stiffness Matrix

% Compliance matrix
matrices.S = [(1/ lamina.E1) (-lamina.v12 / lamina.E1) 0; ...
              (-lamina.v12 / lamina.E1) (1 / lamina.E2) 0; ...
              0 0 (1 / lamina.G12)];

% Stiffness matrix
matrices.Q = inv(matrices.S);

for ii = 1:length(layup)
    
    matrices.Tepsilon(:, :, ii) = [cosd(layup(ii))^2 sind(layup(ii))^2 ...
        ((cosd(layup(ii))) * (sind(layup(ii)))); sind(layup(ii))^2 ...
        cosd(layup(ii))^2 (-(cosd(layup(ii))) * (sind(layup(ii)))); ...
        (-2 * (cosd(layup(ii))) * (sind(layup(ii)))) ...
        (2 * (cosd(layup(ii))) * (sind(layup(ii)))) ...
        ((cosd(layup(ii))^2) - (sind(layup(ii))^2))];
    
    matrices.Qbar(:, :, ii) = transpose(matrices.Tepsilon(:, :, ii)) * ...
        matrices.Q * matrices.Tepsilon(:, :, ii);
    
        matrices.Tsigma(:, :, ii) = [sind(layup(ii))^2 cosd(layup(ii))^2 ...
        (2 *(cosd(layup(ii))) * (sind(layup(ii)))); cosd(layup(ii))^2 ...
        sin(layup(ii))^2 (-2 * (cosd(layup(ii))) * (sind(layup(ii)))); ...
        (-1 * (cosd(layup(ii))) * (sind(layup(ii)))) ...
        (1 * (cosd(layup(ii))) * (sind(layup(ii)))) ...
        ((cosd(layup(ii))^2) - (sind(layup(ii))^2))];
    
end

%% ABD Matrix

% Initializing ABD matrices
[matrices.A, matrices.B, matrices.D] = deal(zeros(3, 3));

% Solving for ABD matrices
for ii = 1:3
    for jj = 1:3
        for kk = 2:length(layup) + 1
            
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

% ABD Matrix
matrices.ABD = [matrices.A matrices.B; matrices.B matrices.D];

%% Global Hygrothermal Properties & Free Strain
for ii = 1:length(layup)
    
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
for ii = 1:length(layup)

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

%% Midplane deformation
hygrothermal.midplaneDeformation = inv(matrices.ABD) * ...
    [hygrothermal.N; hygrothermal.M];

for ii = 1:length(layup)
    
    hygrothermal.globStrain(:, :, ii) = ...
        hygrothermal.midplaneDeformation(1:3) + ...
        (laminate.zCoord(ii) * hygrothermal.midplaneDeformation(4:end)) - ...
        hygrothermal.freeStrain(:, :, ii);
end

% Hygrothermal stress (lamina) and local strain
for ii = 1:length(layup)

    % Global Stress
    hygrothermal.globStress(:, :, ii) = matrices.Qbar(:, :, ii) * ...
        (hygrothermal.midplaneDeformation(1:3) + (zMod(ii) * ...
        hygrothermal.midplaneDeformation(4:end)) - ...
        hygrothermal.freeStrain(:, :, ii));
    
    % Local Stress
    hygrothermal.locStress(:, :, ii) = matrices.Tsigma(:, :, ii) * ...
        hygrothermal.globStress(:, :, ii);
    
    % Local strain
    hygrothermal.locStrain(:, :, ii) = matrices.Tepsilon(:, :, ii) * ...
        hygrothermal.globStrain(:, :, ii);
    
end
