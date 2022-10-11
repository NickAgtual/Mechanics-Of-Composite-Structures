clear; clc; close all

theta = linspace(0, 90, 270);
%% Material Properties

% All material properties abide by the following order
% [Carbon/Epoxy, Boron/Aluminum, Glass/Epoxy]

% Longitudinal modulus
prop.E1 = [147 235 41] * 10 ^ 3;
% Transverse in-plane modulus
prop.E2 = [10.3 137 10.4] * 10 ^ 3;
% In-plane Poisson's Ratio
prop.v12 = [.27 .3 .28];
% In-plane shear modulus
prop.G12 = [7 47 4.3] * 10 ^ 3;


%% Engineering Moduli

youngsModulus = zeros(1, length(theta), length(prop.E1));
shearModulus = zeros(1, length(theta), length(prop.E1));

for ii = 1:length(prop.E1)
    
    moduli.Ex(:, :, ii) = 1 ./ ((cosd(theta).^2 ./ prop.E1(ii)) .* ...
                (cosd(theta).^2 - (sind(theta).^2 .* prop.v12(ii))) + ...
                (sind(theta).^2 ./ prop.E2(ii)) .* ...
                (sind(theta).^2 - (cosd(theta).^2 .* prop.v12(ii))) + ...
                ((cosd(theta).^2 .* sind(theta).^2) ./ prop.G12(ii)));
            
    youngsModulus(:, :, ii) = moduli.Ex(:, :, ii) ./ prop.E2(ii);

    moduli.Gxy(:, :, ii) = 1 ./ ((((4 * cosd(theta).^2 .* sind(theta).^2) ...
        ./ prop.E1(ii)) .* (1 + prop.v12(ii))) + ...
        (((4 * cosd(theta).^2 .* sind(theta).^2) ...
        ./ prop.E2(ii)) .* (1 + prop.v12(ii))) + ((((cosd(theta).^2) - ...
        (sind(theta).^2)).^2) ./ prop.G12(ii)));
               
   shearModulus(:, :, ii) = 10 * moduli.Gxy(:, :, ii) ./ prop.G12(ii);
   
   moduli.Etasx(:, :, ii) = ((((2 .* cosd(theta).^3 .* sind(theta)) ./ ...
       prop.E1(ii)) .* (1 + prop.v12(ii))) - ...
       (((2 .* cosd(theta) .* sind(theta).^3) ./ prop.E2(ii)) .* ...
       (1 + prop.v12(ii))) - (((cosd(theta) .* sind(theta)) .* ...
       ((cosd(theta).^2 - sind(theta).^2))) ./ prop.G12(ii))) .* ...
       moduli.Gxy(:, :, ii);
   
   moduli.poisson(:, :, ii) = ((((cosd(theta).^2) ./ prop.E1(ii)) .* ...
       (((cosd(theta).^2) .* prop.v12(ii)) - sind(theta).^2)) + ...
       (((sind(theta).^2) ./ prop.E2(ii)) .* ...
       (((sind(theta).^2) .* prop.v12(ii)) - cosd(theta).^2)) + ...
       (((cosd(theta).^2) .* sind(theta).^2) ./ prop.G12(ii))) .* ...
       moduli.Ex(:, :, ii);
                   
end

figure(1)

color = ['r', 'b', 'k'];

displayNamesY = {'Carbon/Epoxy Young''s Modulus' , ...
                'Boron/Aluminum Young''s Modulus', ...
                'Glass/Epozy Young''s Modulus'};
            
displayNamesS = {'Carbon/Epoxy Shear Modulus' , ...
                'Boron/Aluminum Shear Modulus', ...
                'Glass/Epoxy Shear Modulus'};
            
for ii = 1:length(prop.E1)
    
    plot(theta, youngsModulus(:, :, ii), strcat(color(ii), '-.'), ...
        'DisplayName', displayNamesY{ii})
    hold on
    grid on
    grid minor
    plot(theta, shearModulus(:, :, ii), strcat(color(ii), '-'), ...
        'DisplayName', displayNamesS{ii})
end

title('\emph{E $\&$ G as a Function of Fiber Orientation}', 'fontsize', ...
    16, 'Interpreter', 'Latex')
xlabel('\emph{$\theta$ (degress)}', 'fontsize', 14, 'Interpreter', 'Latex')
ylabel('\emph{${\frac{E_x}{E}}$, 10G${_{xy}}$}/{G$_{12}$}',...
    'fontsize', 14, 'Interpreter', 'Latex')
legend('location', 'southoutside', 'NumColumns', 3)

figure(2)

displayNamesE = {'Carbon/Epoxy Shear Coupling Coeff.' , ...
                'Boron/Aluminum Shear Coupling Coeff.', ...
                'Glass/Epoxy Shear Coupling Coeff.'};
            
displayNamesV = {'Carbon/Epoxy Poisson''s Ratio' , ...
                'Boron/Aluminum Poisson''s Ratio', ...
                'Glass/Epoxy Poisson''s Ratio'};

for ii = 1:length(prop.E1)
    
    plot(theta, -moduli.Etasx(:, :, ii), strcat(color(ii), '-.'), ...
        'DisplayName', displayNamesE{ii})
    hold on 
    grid on
    grid minor
    plot(theta, moduli.poisson(:, :, ii), strcat(color(ii), '-'), ...
        'DisplayName', displayNamesV{ii})
    
end

title('\emph{$\eta_{sx}$ $\&$ $\nu_{xy}$ as a Function of Fiber Orientation}', ...
    'fontsize', 16, 'Interpreter', 'Latex')
xlabel('\emph{$\theta$ (degress)}', 'fontsize', 14, 'Interpreter', 'Latex')
ylabel('\emph{-$\eta_{sx}$, $\nu_{xy}$}',...
    'fontsize', 14, 'Interpreter', 'Latex')
legend('location', 'southoutside', 'NumColumns', 3)

