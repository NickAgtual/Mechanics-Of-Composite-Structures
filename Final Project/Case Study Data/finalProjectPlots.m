%% Final Project Plots
clear; clc; close all

%% Reading Excel Data

% File names for each load case
fileNames = {'case1.xlsx', 'case2.xlsx', 'case3.xlsx'};
stressStrainType = {'Longitudinal', 'Transverse', 'Shear'};

% Initializing data cell array
data = cell(1, 3);

for ii = 1:length(fileNames)
    
    % Reading each excel file and storing in 'data' cell array
    data{ii} = readcell(fileNames{ii});
    
end

%% Plots

% ----- Global Stress -----

% Creating new figure
figure(1)

for ii = 1:3
    for jj = 1:length(fileNames)
    
    % Creating subplot
    subplot(3, 1, jj)
    
    plot(cell2mat(data{1, ii}(13 + jj, 2:end)), ...
        cell2mat(data{1, ii}(6, 2:end)), 'DisplayName', ...
        strcat('Case ', num2str(ii)));

    % Plot parameters
    hold on
    grid on
    grid minor
    
    titleText = strcat(stressStrainType{jj}, ' Stress');
    
    % Plot descriptors
    xlabel('\emph {Global Stress (psi)}','fontsize',12, ...
        'Interpreter', 'latex');
    ylabel('\emph {z-Location (in)}','fontsize',12,'Interpreter', ...
        'latex');
    title(titleText,'fontsize',14,'Interpreter',...
        'latex')
    legend('location', 'Best', 'Interpreter', 'latex')
    
    end
end

% Subplot title
sgtitle('\emph {Global Stress at Each Ply}', 'fontsize', 16, ...
    'Interpreter', 'latex')
    

% ----- Global Strain -----

% Creating new figure
figure(2)

for ii = 1:3
    for jj = 1:length(fileNames)
    
    % Creating subplot
    subplot(3, 1, jj)
    
    plot(cell2mat(data{1, ii}(18 + jj, 2:end)), ...
        cell2mat(data{1, ii}(6, 2:end)), 'DisplayName', ...
        strcat('Case ', num2str(ii)));

    % Plot parameters
    hold on
    grid on
    grid minor
    
    titleText = strcat(stressStrainType{jj}, ' Strain');
    
    % Plot descriptors
    xlabel('\emph {Global Strain (psi)}','fontsize',12, ...
        'Interpreter', 'latex');
    ylabel('\emph {z-Location (in)}','fontsize',12,'Interpreter', ...
        'latex');
    title(titleText,'fontsize',14,'Interpreter',...
        'latex')
    legend('location', 'Best', 'Interpreter', 'latex')
    
    end
end

% Subplot title
sgtitle('\emph {Global Strain at Each Ply}', 'fontsize', 16, ...
    'Interpreter', 'latex')