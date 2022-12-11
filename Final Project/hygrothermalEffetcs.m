function [hygrothermal] = ...
    hygrothermalEffetcs(ssMod, hygrothermal, z, Tepsilon, Tsigma, ABD, ...
    Qbar, ss)

% Temperature difference
hygrothermal.deltaT = hygrothermal.Tf - hygrothermal.T0;

%% Global Hygrothermal Properties & Free Strain
for ii = 1:length(ssMod)
    
    % Global hygrothermal properties
    hygrothermal.alphaGlobal(:, :, ii) = inv(Tepsilon(:, :, ii)) * ...
        hygrothermal.alpha';
    
    % Free strain
    hygrothermal.freeStrain(:, :, ii) = ...
        (hygrothermal.alphaGlobal(:, :, ii) .* hygrothermal.deltaT);
    
end

%% Hygrothermal Loads

% Initializing hygrothermal loads
[hygrothermal.Nprelim, hygrothermal.Mprelim] = deal(zeros(3, 1));

% Calculating forces and moments due to hygrothermal conditions
for ii = 1:length(ssMod)
    
    if ii < (length(ss)/2) + 1
        
        % Forces (does not include mutliplication of temp diff.)
        hygrothermal.Nprelim = hygrothermal.Nprelim + ((Qbar(:, :, ii) ...
            * hygrothermal.alphaGlobal(:, :, ii) .* (z(ii+1) - z(ii))));
        
        % Moments (does not include mutliplication of temp diff.)
        hygrothermal.Mprelim = hygrothermal.Mprelim + ((Qbar(:, :, ii) * ...
            hygrothermal.alphaGlobal(:, :, ii) * ((z(ii + 1) ^ 2) ...
            - (z(ii) ^ 2))));
        
    elseif ii == (length(ss)/2) +1
        
        continue
        
    elseif ii > (length(ss)/2) + 1
        
        % Forces (does not include mutliplication of temp diff.)
        hygrothermal.Nprelim = hygrothermal.Nprelim + ((Qbar(:, :, ii) ...
            * hygrothermal.alphaGlobal(:, :, ii) .* (z(ii) - z(ii-1))));
        
        % Moments (does not include mutliplication of temp diff.)
        hygrothermal.Mprelim = hygrothermal.Mprelim + ((Qbar(:, :, ii) * ...
            hygrothermal.alphaGlobal(:, :, ii) * ((z(ii) ^ 2) ...
            - (z(ii-1) ^ 2))));
    end
    
end

% Thermal Loads
hygrothermal.N = hygrothermal.Nprelim * hygrothermal.deltaT;
hygrothermal.M = hygrothermal.Mprelim * hygrothermal.deltaT * .5;

% Hygrothermal midplain deformation
% Since this only accounts for temperature change, the total midplane def. 
% is equivalent to the midplain strain and curverature due to hygrothermal
% conditions
hygrothermal.midplaneDeformation = inv(ABD) * ...
    [hygrothermal.N; hygrothermal.M];

for ii = 1:length(ssMod)
    
    hygrothermal.globStrain(:, :, ii) = ...
        hygrothermal.midplaneDeformation(1:3) + ...
        (z(ii) * hygrothermal.midplaneDeformation(4:end)) - ...
        hygrothermal.freeStrain(:, :, ii);
end

% Hygrothermal stress (lamina) and local strain
for ii = 1:length(ssMod)

    % Global Stress
    hygrothermal.globStress(:, :, ii) = Qbar(:, :, ii) * ...
        (hygrothermal.midplaneDeformation(1:3) + (z(ii) * ...
        hygrothermal.midplaneDeformation(4:end)) - ...
        hygrothermal.freeStrain(:, :, ii));
    
    % Local Stress
    hygrothermal.locStress(:, :, ii) = Tsigma(:, :, ii) * ...
        hygrothermal.globStress(:, :, ii);
    
    % Local strain
    hygrothermal.locStrain(:, :, ii) = Tepsilon(:, :, ii) * ...
        hygrothermal.globStrain(:, :, ii);
    
end

end


