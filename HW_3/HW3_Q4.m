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
[A, B, D] = deal(zeros(3, 3));

for ii = 1:3
    for jj = 1:3
        for kk = 2:length(t)+1
            
            A(ii, jj) = A(ii, jj) + (Qbar(ii, jj, kk-1) * (z(kk) - z(kk-1)));
            
            B(ii, jj) = B(ii, jj) + (.5 * (Qbar(ii, jj, kk-1) * ...
                (z(kk)^2 - z(kk-1)^2)));
            
            D(ii, jj) = D(ii, jj) + ((1/3) * (Qbar(ii, jj, kk-1) * ...
                (z(kk)^3 - z(kk-1)^3)));
            
        end
    end
end
    