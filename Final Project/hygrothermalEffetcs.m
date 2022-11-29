function [hygrothermal] = ...
    hygrothermalEffetcs(ss, hygrothermal, z, Te, ABD, Qbar)

% Temperature difference
hygrothermal.deltaT = hygrothermal.Tf - hygrothermal.T0;

% Moisture difference
hygrothermal.deltaC = hygrothermal.Cf - hygrothermal.C0;


%% Local Hygrothermal Properties & Free Strain
for ii = 1:length(ss)
    
    % Local hygrothermal properties
    hygrothermal.alphaLocal(:, :, ii) = Te(:, :, ii) * hygrothermal.alpha';
    hygrothermal.betaLocal(:, :, ii) = Te(:, :, ii) * hygrothermal.beta';
    
    % Free strain
    hygrothermal.freeStrain(:, :, ii) = ...
        (hygrothermal.alphaLocal(:, :, ii) .* hygrothermal.deltaT) + ...
        (hygrothermal.betaLocal(:, :, ii) .* hygrothermal.deltaC);
    
end

%% Hygrothermal Loads

% Initializing hygrothermal loads
[hygrothermal.N, hygrothermal.M] = deal(zeros(3, 1));

for ii = 2:length(ss) + 1

    hygrothermal.N = (hygrothermal.N + (Qbar(:, :, ii - 1) ...
        * hygrothermal.alphaLocal(:, :, ii - 1) .* (z(ii) - z(ii - 1)))); 
    % ACCOUNT FOR MOISTURE TOO
    
    hygrothermal.M = (hygrothermal.M + (Qbar(:, :, ii - 1) ...
        * hygrothermal.betaLocal(:, :, ii - 1) .* ((z(ii) ^ 2) - ...
        (z(ii - 1) ^ 2)))); % ACCOUNT FOR MOISTURE TOO
end

end


