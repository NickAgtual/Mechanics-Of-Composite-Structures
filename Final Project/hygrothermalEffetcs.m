function [hygrothermal] = ...
    hygrothermalEffetcs(ss, hygrothermal, z, Te, Tsigma, ABD, Qbar)

% Temperature difference
hygrothermal.deltaT = hygrothermal.Tf - hygrothermal.T0;

%% Global Hygrothermal Properties & Free Strain
for ii = 1:length(ss)
    
    % Global hygrothermal properties
    hygrothermal.alphaGlobal(:, :, ii) = Te(:, :, ii) * hygrothermal.alpha';
    
    % Free strain
    hygrothermal.freeStrain(:, :, ii) = ...
        (hygrothermal.alphaGlobal(:, :, ii) .* hygrothermal.deltaT);
    
end

%% Hygrothermal Loads

% Initializing hygrothermal loads
[hygrothermal.N, hygrothermal.M] = deal(zeros(3, 1));

% Calculating forces and moments due to hygrothermal conditions
for ii = 1:length(ss)

    % Forces (does not include mutliplication of temp diff.)
    hygrothermal.Nprelim = (hygrothermal.N + (Qbar(:, :, ii) ...
        * hygrothermal.alphaGlobal(:, :, ii) .* (z(ii+1) - z(ii)))); 
    
    % Moments (does not include mutliplication of temp diff.)
    hygrothermal.Mprelim = (hygrothermal.M + (Qbar(:, :, ii) * ...
        hygrothermal.alphaGlobal(:, :, ii) * ((z(ii + 1) ^ 2) ...
        - (z(ii) ^ 2))));
end

% Thermal Loads
hygrothermal.N = hygrothermal.Nprelim * hygrothermal.deltaT;
hygrothermal.M = hygrothermal.Mprelim * hygrothermal.deltaT * .5;

% Hygrothermal global strain (lamina)
% Since this only accounts for temperature change, the global strain here
% is equivalent to the midplain strain and curverature due to hygrothermal
% conditions
hygrothermal.globStrain = inv(ABD) * [hygrothermal.N; hygrothermal.M];

% Hygrothermal stress (lamina)
for ii = 1:length(ss)

    % Global Stress
    hygrothermal.globStress(:, :, ii) = Qbar(:, :, ii) * ...
        (hygrothermal.globStrain(1:3) + (z(ii) * ...
        hygrothermal.globStrain(4:end)) - ...
        hygrothermal.freeStrain(:, :, ii));
    
    % Local Stress
    hygrothermal.localStress(:, :, ii) = Tsigma(:, :, ii) * ...
        hygrothermal.globStress(:, :, ii);
end

end


