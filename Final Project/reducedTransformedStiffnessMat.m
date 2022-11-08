function [Qbar, S, Te] = reducedTransformedStiffnessMat(moduli, ss)

%% Transformed Reduced Stiffness Matrix

% Reduced compliance matrix

S = [(1/moduli.E1) (-moduli.v12/moduli.E1) 0;
    (-moduli.v12/moduli.E1) (1/moduli.E2) 0;
    0 0 (1/moduli.G12)];

% Reduced stiffness matrix
Q = inv(S);

% Reduced strain transformation matrix

% Initializing Te and Qbar matrix
Te = zeros(3, 3, length(ss));
Qbar = zeros(3, 3, length(ss));

for ii = 1:length(ss)
    
    % Reduced strain transformation matrix
    Te(:,:, ii) = [(cosd(ss(ii)).^2) (sind(ss(ii)).^2) (cosd(ss(ii)) .* sind(ss(ii)));
        (sind(ss(ii)).^2) (cosd(ss(ii)).^2) -(cosd(ss(ii)) .* sind(ss(ii)));
        (-2 .* (cosd(ss(ii)) .* sind(ss(ii)))) ...
        (2 .* (cosd(ss(ii)) .* sind(ss(ii)))) ...
        ((cosd(ss(ii)).^2) - (sind(ss(ii)).^2))];
    
    Qbar(:, :, ii) = transpose(Te(:, :, ii)) * Q * Te(:, :, ii);
    
end

end