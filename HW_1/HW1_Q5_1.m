%% Defining Global System

globSys.x = [1; 0; 0];
globSys.y = [0; 1; 0];
globSys.z = [0; 0; 0];

% State of stress [sigmaX, sigmaY, sigmaXY]
globSys.stress = [80 -50 0 0 0 15]; % ksi

%% Defining Local System

theta = 20;

localSys.x = [cosd(theta); sind(theta); 0];
localSys.y = [-sind(theta); cosd(theta); 0];
localSys.z = [0; 0; 1];

%% Plane Stress and Orthotropic Transformation

reducedStressTransformation = ...
    [cosd(theta)^2 sind(theta)^2 (2 * cosd(theta) * sind(theta));
    sind(theta)^2 cosd(theta)^2 (-2 * cosd(theta) * sind(theta));
    (-cosd(theta) * sind(theta)) (cosd(theta) * sind(theta)) ...
    (cosd(theta)^2 - sind(theta)^2)];

localStress = reducedStressTransformation * [80 -50 15]';


