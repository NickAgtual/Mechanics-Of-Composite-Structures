%% HW3 Q4
clear; clc; close all

%% Material Properties

% Longitudinal modulus
lamina.E1 = 42.7 * 10 ^ 6; % psi

% Transverse moudulus
lamina.E2 = .92 * 10 ^ 6; % psi

% In-plane shear modulus
lamina.G12 = .71 * 10 ^ 6; % psi

% Poisson's ratio
lamina.v12 = .23;

% Stacking sequences
laminate.layup = [45 90 30 0];

% Thickness of lamina
lamina.thickness = .0008; % in

% Lamina coordinates within laminate

z = zeros(1, length(laminate.layup) / 2);

for ii = 1: length(laminate.layup) / 2
    
    z(ii) = lamina.thickness * ii;
    
end

laminate.zCoord = [-z 0 z];

%% 1 Calculate the A, B, D, a, b, d Matrices for the Layup [45/90/30/0]

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
        cosd(laminate.layup(ii))^2 - sind(laminate.layup(ii))];
    
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

% a matrix
matrices.a = matrices.abd(1:3, 1:3);

% b matrix
matrices.b = matrices.abd(1:3, 4:end);

% c matrix
matrices.d = matrices.abd(4:end, 4:end);

    