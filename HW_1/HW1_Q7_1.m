%% Material Properties

% Longitudinal modulus
prop.E1 = 140 * 10^3; % MPa
% Transverse in-plane modulus
prop.E2 = 10 * 10^3; % MPa

prop.G12 = 7 * 10^3; % MPa

prop.v12 = .3;

%% Global Stress State

syms tau strainXX strainYY

globSys.stress = [10 0 tau]; %MPa

globSys.strain = [strainXX strainYY 0];

%% Local System

theta = 45;

reducedStressTransformation = ...
    [cosd(theta)^2 sind(theta)^2 (2 * cosd(theta) * sind(theta));
    sind(theta)^2 cosd(theta)^2 (-2 * cosd(theta) * sind(theta));
    (-cosd(theta) * sind(theta)) (cosd(theta) * sind(theta)) ...
    (cosd(theta)^2 - sind(theta)^2)];

complianceMat = [(1 / prop.E1) (-prop.v12 / prop.E1) 0;
                  (-prop.v12 / prop.E1) (1 / prop.E2) 0;
                  0 0 (1 / prop.G12)];
              
% Transformed reduced compliance matrix 
sBar = transpose(reducedStressTransformation) * complianceMat * ...
    reducedStressTransformation;

globSys.strain = sBar * globSys.stress';

tau = double(solve(globSys.strain(3) == 0, tau));

%% ---------------------------------- ALTERNATE METHOD -------------------
% global stress --> local Stress --> local strain --> global strain  


