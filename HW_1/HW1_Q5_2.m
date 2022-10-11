%% Defining Global System

globSys.x = [1; 0; 0];
globSys.y = [0; 1; 0];
globSys.z = [0; 0; 0];

%% Defining Local System

theta = -30;

localSys.x = [cosd(theta); sind(theta); 0];
localSys.y = [-sind(theta); cosd(theta); 0];
localSys.z = [0; 0; 1];

% Plane Stress [sigmaX sigmaY sigmaXY]'
localSys.stress = [2280 0 0]'; % MPa

%% Material Properties

% Longitudinal modulus
prop.E1 = 147 * 10^3; % MPa
% Transverse in-plane modulus
prop.E2 = 10.3 * 10^3; %MPa

% Major in-palne poisson's ratio
prop.v12 = .27;
% Out-of-plane poisson's ratio
prop.v23 = .54;
% In-plane shear modulus
prop.G12 = 7 * 10^3; % MPa


%% Stress Transformation
% Going from local to global therefore * -1
theta = -30 * -1;

reducedStressTransformation = ...  
    [cosd(theta)^2 sind(theta)^2 (2 * cosd(theta) * sind(theta));
    sind(theta)^2 cosd(theta)^2 (-2 * cosd(theta) * sind(theta));
    (-cosd(theta) * sind(theta)) (cosd(theta) * sind(theta)) ...
    (cosd(theta)^2 - sind(theta)^2)];

globSys.stress = reducedStressTransformation * localSys.stress;

