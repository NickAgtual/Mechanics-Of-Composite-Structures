function stressStrainPlots(superimposedParam, z)

stressStrainType = {'Longitudinal', 'Transverse', 'Shear'};

% Creating vector with magnitude of stress and strain
for ii = 1:length(stressStrainType)
    for jj = 1:length(z)
    
    % Concatenating local stress values (3 filed for each stress type)
    localStress.(stressStrainType{ii})(jj) = ...
        superimposedParam.locStress(ii, 1, jj);
    
    % Concatenating local strain values (3 fields for each strain type)
    localStrain.(stressStrainType{ii})(jj) = ...
        superimposedParam.locStrain(ii, 1, jj);
    end
end

%% Local Stress Plot
for ii = 1:length(stressStrainType)
    
% Creating new figure
figure(ii)

% Plotting stress vs. z-coordinate
plot(localStress.(stressStrainType{ii}), z, '-o')

% Plot parameters
grid on
grid minor

% Plot descriptors
xlabel('\emph {Global Stress ()}','fontsize',14,'Interpreter',...
    'latex');
ylabel('\emph {z-Location (in)}','fontsize',14,'Interpreter','latex');
title('\emph {Global Stress in Laminate}','fontsize',16,'Interpreter',...
    'latex')
legend('location', 'Best', 'Interpreter', 'latex')

end


%% Local Strain Plot
for ii = 1:length(stressStrainType)
% Creating new figure
figure(ii + 3)

% Plotting strain vs. z-coordinate
plot(localStrain.(stressStrainType{ii}), z, '-o')

% Plot parameters
grid on
grid minor

% Plot descriptors
xlabel('\emph {Global Strain ()}','fontsize',14,'Interpreter',...
    'latex');
ylabel('\emph {z-Location (in)}','fontsize',14,'Interpreter','latex');
title('\emph {Global Strain in Lamina}','fontsize',16,'Interpreter',...
    'latex')
legend('location', 'Best', 'Interpreter', 'latex')

end
end
