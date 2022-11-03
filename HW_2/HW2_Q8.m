%% Material Properties

% Longitudinal modulus
prop.E1 = 140; % GPa

% Transverse modulus
prop.E2 = 10; % GPa

% Poisson's Ratio
prop.v12 = .3;

% Shear modulus
prop.G12 = 7; % GPa

% Ply thickness
prop.t = .127 * 10 ^ -3; % m

% Layup sequence
t = [0 0 90 90];

z = [-.254 -.127 0 .127 .254] * 10 ^ -3;

% Loading
Nx = 5000;

%% Transformed Reduced Stiffness Matrix

% Reduced compliance matrix

S = [(1/prop.E1) (-prop.v12/prop.E1) 0;
    (-prop.v12/prop.E1) (1/prop.E2) 0;
    0 0 (1/prop.G12)];

% Reduced stiffness matrix
Q = inv(S);

% Reduced strain transformation matrix

% Initializing Te and Qbar matrix
Te = zeros(3, 3, length(t));
Qbar = zeros(3, 3, length(t));

for ii = 1:length(t)
    
    % Reduced strain transformation matrix
    Te(:,:, ii) = [(cosd(t(ii)).^2) (sind(t(ii)).^2) (cosd(t(ii)) .* sind(t(ii)));
        (sind(t(ii)).^2) (cosd(t(ii)).^2) -(cosd(t(ii)) .* sind(t(ii)));
        (-2 .* (cosd(t(ii)) .* sind(t(ii)))) ...
        (2 .* (cosd(t(ii)) .* sind(t(ii)))) ...
        ((cosd(t(ii)).^2) - (sind(t(ii)).^2))];
    
    Qbar(:, :, ii) = transpose(Te(:, :, ii)) * Q * Te(:, :, ii);
    
end

%% ABD Matrix

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


%% Deflection caused by Uniaxial Load

% Laminate stiffness matrix
ABD = [A B;B D];

% Uniaxial loading
uniaxialLoading = [Nx 0 0]';

% Deflection due to uniaxial loading
strain.uniaxial = inv(A) * uniaxialLoading;

%% Solving for Moments

moments = B * strain.uniaxial;

