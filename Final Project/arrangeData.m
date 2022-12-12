function [toExport] = arrangeData(superimposedParam, z, ...
    deformationAtMidplain, hygrothermal, loading, inputStruct, GUI)

%% Rearranging Existing Data

% Reshaping local stress & strain vectors to 3 x n matrices
toExport.locStress = reshape(superimposedParam.locStress, [3, length(z)]);
toExport.locStrain = reshape(superimposedParam.locStrain, [3, length(z)]);

% Separating midplane deformation
toExport.midplaneStrain = deformationAtMidplain(1:3);
toExport.midplaneCurverature = deformationAtMidplain(4:6);

% Total Loads
toExport.totalForces = hygrothermal.N + loading(1:3)';
toExport.totalMoments = hygrothermal.M + loading(4:6)';

%% Creating Cell Array to Write to Excel

% Writing total forces and moments
for ii = 1:6
    
    % Writing forces
    if (1 <= ii) && (ii <= 3)
        
        toExport.write{ii, 1} = strcat(inputStruct.loading{ii}, 'Tot');
        toExport.write{ii, 2} = toExport.totalForces(ii);
        
    % Writing moments
    elseif (4 <= ii) && (ii <= 6)
        
        toExport.write{ii - 3,4} = strcat(inputStruct.loading{ii}, 'Tot');
        toExport.write{ii - 3, 5} = toExport.totalMoments(ii - 3);
        
    end
end

% Writing change in temperature
toExport.write{1, 7} = 'DeltaT';
toExport.write{1, 8} = hygrothermal.deltaT';

% Writing ply number
toExport.write{5, 1} = 'Ply';
for ii = 1:length(z)
    toExport.write{5, 1 + ii} = ii;
end

% Writing z-coordinate
toExport.write{6, 1} = 'z';
for ii = 1:length(z)
    toExport.write{6, 1 + ii} = z(ii);
end

% Writing midplane deformations
toExport.write{8, 1} = 'Midplane Deformation';
midplaneNames = {'eps0x', 'eps0y', 'eps0xy', 'kx', 'ky', 'kxy'};

for ii = 1:6
    
    % Writing midplain strain
    if (1 <= ii) && (ii <= 3)
        
        toExport.write{ii + 8, 1} = midplaneNames{ii};
        toExport.write{ii + 8, 2} = toExport.midplaneStrain(ii);
        
    % Writing midplane curveratures
    elseif (4 <= ii) && (ii <= 6)
        
        toExport.write{ii + 5, 4} = midplaneNames{ii};
        toExport.write{ii + 5, 5} = toExport.midplaneCurverature(ii - 3);
        
    end
end

% Writing local stress
toExport.write{13, 1} = 'Local Stress';
stressNames = {'sigma1', 'sigma2', 'sigma12'};

for ii = 1:length(stressNames)
    for jj = 1:length(toExport.locStress) + 1
        
        if jj == 1
            
            toExport.write{ii + 13, 1} = stressNames{ii};
            
        elseif jj > 1
            
            toExport.write{ii + 13, jj} = toExport.locStress(ii, jj - 1);
            
        end
    end
    
end

% Writing local strain
toExport.write{18, 1} = 'Local Strain';
strainNames = {'epsilon1', 'epsilon2', 'gamma12'};

for ii = 1:length(strainNames)
    for jj = 1:length(toExport.locStress) + 1
        
        if jj == 1
            
            toExport.write{ii + 18, 1} = strainNames{ii};
            
        elseif jj > 1
            
            toExport.write{ii + 18, jj } = toExport.locStrain(ii, jj - 1);
            
        end
    end
    
end

%% Exporting Data to Excel Sheet

if GUI == 0
    
    % Question dialog box
    answer = questdlg('Export Data to Excel?', 'Export', 'Yes', 'No', 'Yes');
    
    % Handle response
    switch answer
        case 'Yes'
            % Creating file name using current date and time
            dateTime = datestr(now, 'ddmmyy-HHMM');
            toExport.fileName = strcat('Composite Analysis - ', dateTime, ...
                '.xlsx');
            
            % Writing to excel file
            writecell(toExport.write,toExport.fileName)
            
        case 'No'
            delete(answer)
            
    end
    
end

end




