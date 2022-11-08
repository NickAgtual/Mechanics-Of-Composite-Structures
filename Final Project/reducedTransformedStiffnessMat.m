function [Qbar] = reducedTransformedStiffnessMat(moduli)

%% Transformed Reduced Stiffness Matrix

% Reduced compliance matrix

S = [(1/moduli.E1) (-moduli.v12/moduli.E1) 0;
    (-moduli.v12/moduli.E1) (1/moduli.E2) 0;
    0 0 (1/moduli.G12)];

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

end