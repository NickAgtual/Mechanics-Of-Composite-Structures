function [x] = stressStrainPlots(superimposedParam, zMod)

% Creating vector with magnitude of stress and strain
for ii = 1:length(superimposedParam.globStress)
    globStressVec(ii) = norm(superimposedParam.globStress(:, ii));
    globStrainVec(ii) = norm(superimposedParam.globStrain(:, ii));
end

%% Global Stress Plot

% Creating new figure
figure(1)

% Plotting stress vs. z-coordinate
plot(globStressVec, zMod, '-o')

grid on
grid minor

% Plot descriptors
xlabel('\emph {Global Stress ()}','fontsize',14,'Interpreter',...
    'latex');
ylabel('\emph {z-Location (in)}','fontsize',14,'Interpreter','latex');
title('\emph {Global Stress in Laminate}','fontsize',16,'Interpreter',...
    'latex')
legend('location', 'Best', 'Interpreter', 'latex')

%% Global Strain Plot

% Creating new figure
figure(2)

% Plotting strain vs. z-coordinate
plot(globStrainVec, zMod, '-o')

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
