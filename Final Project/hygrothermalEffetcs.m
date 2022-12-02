function [hygrothermal] = ...
    hygrothermalEffetcs(ss, hygrothermal, z, Te, ABD, Qbar)

% Temperature difference
hygrothermal.deltaT = hygrothermal.Tf - hygrothermal.T0;

%% Global Hygrothermal Properties & Free Strain
for ii = 1:length(ss)
    
    % Global hygrothermal properties
    hygrothermal.alphaLocal(:, :, ii) = Te(:, :, ii) * hygrothermal.alpha';
    
    % Free strain
    hygrothermal.freeStrain(:, :, ii) = ...
        (hygrothermal.alphaLocal(:, :, ii) .* hygrothermal.deltaT);
    
end

%% Hygrothermal Loads

% Initializing hygrothermal loads
[hygrothermal.N, hygrothermal.M] = deal(zeros(3, 1));

for ii = 2:length(ss) + 1

    hygrothermal.N = (hygrothermal.N + (Qbar(:, :, ii - 1) ...
        * hygrothermal.alphaLocal(:, :, ii - 1) .* (z(ii) - z(ii - 1)))); 
    
    hygrothermal.M = (hygrothermal.M + (Qbar(:, :, ii - 1)));
end

N = hygrothermal.N
M = hygrothermal.M
% Hygrothermal global Strain



% Hygrothermal global stress 

end


