function [Sbar, Tsigma] = reducedTransformedComplianceMat(S, ssMod)

% Initializing Te and Qbar matrix
[Tsigma, Sbar] = deal(zeros(3, 3, length(ssMod)));

for ii = 1:length(ssMod)
    
    % Reduced strain transformation matrix
    Tsigma(:,:, ii) = [(cosd(ssMod(ii)).^2) (sind(ssMod(ii)).^2) ...
        (2 * cosd(ssMod(ii)) .* sind(ssMod(ii)));
        (sind(ssMod(ii)).^2) (cosd(ssMod(ii)).^2) ((-2 * cosd(ssMod(ii)) ...
        .* sind(ssMod(ii)))); (-1 * (cosd(ssMod(ii)) .* sind(ssMod(ii)))) ...
        ((cosd(ssMod(ii)) .* sind(ssMod(ii)))) ...
        ((cosd(ssMod(ii)).^2) - (sind(ssMod(ii)).^2))];
    
    % Reduced transformed compliance matrix
    Sbar(:, :, ii) = transpose(Tsigma(:, :, ii)) * S * Tsigma(:, :, ii);
    
end

end
