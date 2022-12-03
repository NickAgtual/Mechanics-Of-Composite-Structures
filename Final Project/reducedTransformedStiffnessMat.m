function [Qbar, S, Tepsilon] = reducedTransformedStiffnessMat(moduli, ss)

%% Transformed Reduced Stiffness Matrix

% Reduced compliance matrix
S = [(1/moduli.E1) (-moduli.v12/moduli.E1) 0;
    (-moduli.v12/moduli.E1) (1/moduli.E2) 0;
    0 0 (1/moduli.G12)];

% Reduced stiffness matrix
Q = inv(S);

% Initializing Te and Qbar matrix
[Tepsilon, Qbar] = deal(zeros(3, 3, length(ss)));

for ii = 1:length(ss)
    
    % Reduced strain transformation matrix
    Tepsilon(:,:, ii) = [(cosd(ss(ii)).^2) (sind(ss(ii)).^2) ...
        (cosd(ss(ii)) .* sind(ss(ii)));
        (sind(ss(ii)).^2) (cosd(ss(ii)).^2) -(cosd(ss(ii)) ...
        .* sind(ss(ii))); (-2 .* (cosd(ss(ii)) .* sind(ss(ii)))) ...
        (2 .* (cosd(ss(ii)) .* sind(ss(ii)))) ...
        ((cosd(ss(ii)).^2) - (sind(ss(ii)).^2))];
    
    % Reduced trainsformed stiffness matrix
    Qbar(:, :, ii) = transpose(Tepsilon(:, :, ii)) * Q * ...
        Tepsilon(:, :, ii);
    
end

end