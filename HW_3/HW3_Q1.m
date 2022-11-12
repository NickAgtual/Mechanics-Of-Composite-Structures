%% HW3 Q1
clear; clc; close all;

%% Material Properties

% Longitudinal modulus
prop.E1 = 147; % GPa

% Transverse modulus
prop.E2 = 10.3; % GPa

% Shear modulus
prop.G12 = 7; % GPa

% Poisson's ratio
prop.v12 = .27;

layup = [30 30 0 30 30];

%% Calculating [a] Matrix

% Compliance matrix
S = [(1/ prop.E1) (-prop.v12 / prop.E1) 0;
     (-prop.v12 / prop.E1) (1 / prop.E2) 0;
     0 0 (1 / prop.G12)];
 
% Stiffness matrix
Q = inv(S);

% Initializing Qbar
Qbar = zeros(3, 3, length(layup));

for ii = 1:length(layup)
    
    % Strain transfromation matrix
    strainTransform = [cosd(layup(ii))^2 sind(layup(ii))^2 ...
        ((cosd(layup(ii))^2) * (sind(layup(ii))^2)); sind(layup(ii))^2 ...
        cosd(layup(ii))^2 (-(cosd(layup(ii))^2) * (sind(layup(ii))^2)); ...
        (-2 * (cosd(layup(ii))^2) * (sind(layup(ii))^2)) ...
        (2 * (cosd(layup(ii))^2) * (sind(layup(ii))^2)) ...
        ((cosd(layup(ii))^2) - (sind(layup(ii))^2))];
    
    % Reduced transformed stiffness matrix
    Qbar(:, :, ii) = transpose(strainTransform) * Q * strainTransform;
    
end

% Initializing A matrix
A = zeros(3, 3);

% Solving for A matrix
for ii = 1:3
    for jj = 1:3
        for kk = 2:length(layup)
            
            A(ii, jj) = A(ii, jj) + (Qbar(ii, jj, kk-1) * 1);
            
        end
    end
end

% Solving for a matrix
a = inv(A);

