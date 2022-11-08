function [Sbar] = reducedTransformedComplianceMat(S, ss)


% Initializing Te and Qbar matrix
Tsigma = zeros(3, 3, length(ss));
Sbar = zeros(3, 3, length(ss));

for ii = 1:length(ss)
    
    % Reduced strain transformation matrix
    Tsigma(:,:, ii) = [(cosd(ss(ii)).^2) (sind(ss(ii)).^2) (2 * cosd(ss(ii)) .* sind(ss(ii)));
        (sind(ss(ii)).^2) (cosd(ss(ii)).^2) -(2 * cosd(ss(ii)) .* sind(ss(ii)));
        ((cosd(ss(ii)) .* sind(ss(ii)))) ...
        ((cosd(ss(ii)) .* sind(ss(ii)))) ...
        ((cosd(ss(ii)).^2) - (sind(ss(ii)).^2))];
    
    Sbar(:, :, ii) = transpose(Tsigma(:, :, ii)) * S * Tsigma(:, :, ii);
    
end

end
