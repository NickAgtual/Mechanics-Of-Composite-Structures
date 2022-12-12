function [plotStress, plotStrain] = stressStrainPlots(superimposedParam, ...
    z, GUI)

stressStrainType = {'Longitudinal', 'Transverse', 'Shear'};

% Creating vector with magnitude of stress and strain
for ii = 1:length(stressStrainType)
    for jj = 1:length(z)
        
        % Concatenating local stress values (3 filed for each stress type)
        plotStress.(stressStrainType{ii})(jj) = ...
            superimposedParam.locStress(ii, 1, jj);
        
        % Concatenating local strain values (3 fields for each strain type)
        plotStrain.(stressStrainType{ii})(jj) = ...
            superimposedParam.locStrain(ii, 1, jj);
    end
end

if GUI == 0
    
    %% Local Stress Plot
    for ii = 1:length(stressStrainType)
        
        % Creating new figure
        figure(1)
        
        % Defining subplot position
        subplot(3, 1, ii)
        
        % Plotting stress vs. z-coordinate
        plot(plotStress.(stressStrainType{ii}), z, '-o')
        
        % Plot parameters
        grid on
        grid minor
        
        titleText = strcat(stressStrainType{ii}, ' Stress');
        
        % Plot descriptors
        xlabel('\emph {Global Stress (psi)}','fontsize',12, ...
            Interpreter', 'latex');
        ylabel('\emph {z-Location (in)}','fontsize',12,'Interpreter', ...
            'latex');
        title(titleText,'fontsize',14,'Interpreter',...
            'latex')
        legend('location', 'Best', 'Interpreter', 'latex')
        
    end
    
    % Subplot title
    sgtitle('\emph {Local Stress at Each Ply}', 'fontsize', 16, ...
        'Interpreter', 'latex')
    
    
    %% Local Strain Plot
    for ii = 1:length(stressStrainType)
        % Creating new figure
        figure(2)
        
        % Definign subplot position
        subplot(3, 1, ii)
        
        % Plotting strain vs. z-coordinate
        plot(plotStrain.(stressStrainType{ii}), z, '-o')
        
        % Plot parameters
        grid on
        grid minor
        
        titleText = strcat(stressStrainType{ii}, ' Strain');
        
        % Plot descriptors
        xlabel('\emph {Global Strain}','fontsize',12,'Interpreter',...
            'latex');
        ylabel('\emph {z-Location (in)}','fontsize',12,'Interpreter', ...
            'latex');
        title(titleText,'fontsize',14,'Interpreter',...
            'latex')
        legend('location', 'Best', 'Interpreter', 'latex')
        
    end
    
    % Subplot title
    sgtitle('\emph {Local Strain at Each Ply}', 'fontsize', 16, ...
        'Interpreter', 'latex')
    
end
end
