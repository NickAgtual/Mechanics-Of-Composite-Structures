function [faliure] = faliureCriteria(strength, superimposedParam)

% Listing criterion
criteria = {'Max Stres Criterion', 'Max Strain Criterion', ...
    'Tsai-Hill Criterion', 'Tsai-Wu Criterion'};

% Pupulating cell array with all faliure criterion functions
functions = {@maxStressCriterion, ...
             @maxStrainCriterion, ...
             @tsaiHillCriterion, ...
             @tsaiWuCriterion};

for ii = 1:length(criteria)
    
    % Adding criterion names to faliure structure
    faliure(ii).criterion = criteria{ii};
    
    [faliure(ii).fail, faliure(ii).mode] = ...
        functions{ii}(strength, superimposedParam);
    
end

    function [fail, mode] = maxStressCriterion(strength, ...
            superimposedParam)
        
        % Initializing 'fail' and 'mode' cell arrays
        [fail, mode] = deal(cell(1, length(superimposedParam.locStress)));
        
        % Looping through each ply
        for jj = 1:length(superimposedParam.locStress)
            
            % Checking for fiber faliure
            if (superimposedParam.locStress(1, 1, jj) >= strength.X) || ...
                    (superimposedParam.locStress(2, 1, jj) >= strength.Y)
                
                % Does laminate fail?
                fail{jj} = 'Yes';
                
                % What is the mode of faliure?
                mode{jj} = 'Fiber Faliure';
                
                
                % Checking for matrix faliure (Method 1)
            elseif (superimposedParam.locStress(1, 1, jj) ...
                    <= strength.Xprime) || ...
                    (superimposedParam.locStress(2, 1, jj) <= ...
                    strength.Yprime)
                
                % Does laminate fail?
                fail{jj} = 'Yes';
                
                % What is the mode of faliure?
                mode{jj} = 'Matrix ';
                
            elseif abs(superimposedParam.locStress(3, 1, jj)) >= strength.S
                
                % Does laminate fail?
                fail{jj} = 'Yes';
                
                % What is the mode of faliure?
                mode{jj} = 'Matrix ';
                
            else
                
                % Does laminate fail?
                fail{jj} = 'No';
                
                % What is the mode of faliure?
                mode{jj} = 'N/A';
            end
        end
    end

    function [fail, mode] = maxStrainCriterion(strength, ...
            superimposedParam)
        
                % Initializing 'fail' and 'mode' cell arrays
        [fail, mode] = deal(cell(1, length(superimposedParam.locStress)));
        
                % Looping through each ply
        for jj = 1:length(superimposedParam.locStress)
            
            % Checking for fiber faliure
            if (superimposedParam.locStrain(1, 1, jj) >= ...
                    strength.Xe) || ...
                    (superimposedParam.locStrain(2, 1, jj) >= ...
                    strength.Ye)
                
                % Does laminate fail?
                fail{jj} = 'Yes';
                
                % What is the mode of faliure?
                mode{jj} = 'Fiber Faliure';
                
                % Checking for matrix faliure (Method 1)
            elseif (superimposedParam.locStrain(1, 1, jj) ...
                    <= strength.XePrime) || ...
                    (superimposedParam.locStrain(2, 1, jj) <= ...
                    strength.Yprime)
                
                % Does laminate fail?
                fail{jj} = 'Yes';
                
                % What is the mode of faliure?
                mode{jj} = 'Matrix';
                
            elseif abs(superimposedParam.locStrain(3, 1, jj)) >= ...
                    strength.Se
                
                % Does laminate fail?
                fail{jj} = 'Yes';
                
                % What is the mode of faliure?
                mode{jj} = 'Matrix';
                
            else
                
                % Does laminate fail?
                fail{jj} = 'No';
                
                % What is the mode of faliure?
                mode{jj} = 'N/A';
                
            end
        end
    end

    function [fail, mode] = tsaiHillCriterion(strength, superimposedParam)
        
                % Initializing 'fail' and 'mode' cell arrays
        [fail, mode] = deal(cell(1, length(superimposedParam.locStress)));
        
                % Looping through each ply
        for jj = 1:length(superimposedParam.locStress)
            
            % Variable to compare
            comparisonVar = ...
                ((superimposedParam.locStress(1, 1, jj) / ...
                strength.X) ^ 2) + ...
                ((superimposedParam.locStress(2, 1, jj) / ...
                strength.Y) ^ 2) - ...
                ((superimposedParam.locStress(1, 1, jj) / strength.X) * ...
                (superimposedParam.locStress(2, 1, jj) / strength.X) + ...
                ((superimposedParam.locStress(3, 1, jj) / ...
                strength.S) ^ 2));
            
            if comparisonVar >= 1
                
                % Does laminate fail?
                fail{jj} = 'Yes';
                
            else
                
                % Does laminate fail?
                fail{jj} = 'No';
                
                % What is the mode of faliure?
                mode{jj} = 'N/A';
                
            end    
        end
    end

    function [fail, mode] = tsaiWuCriterion(strength, superimposedParam)
        
                % Initializing 'fail' and 'mode' cell arrays
        [fail, mode] = deal(cell(1, length(superimposedParam.locStress)));
        
                % Looping through each ply
        for jj = 1:length(superimposedParam.locStress)
            
            % Tsai-Wu criterion constants
            constants.F11 = -1 / (strength.Xprime * strength.X);
            constants.F1 = (1 / strength.X) + (1 / strength.Xprime);
            constants.F22 = -1 / (strength.Yprime * strength.Y);
            constants.F2 = (1 / strength.Y) + (1 / strength.Yprime);
            constants.F12 = 1 / (2 * strength.Xprime * strength.X);
            constants.F66 = 1 / (strength.S ^ 2);
            
            % Varibel to compare
            comparisonVar = ...
                (constants.F1 * ...
                superimposedParam.locStress(1, 1, jj)) + ...
                (constants.F2 * ...
                superimposedParam.locStress(2, 1, jj)) + ...
                (constants.F11 * ...
                superimposedParam.locStress(1, 1, jj)^ 2) + ...
                (constants.F22 * ...
                superimposedParam.locStress(2, 1, jj)^ 2) + ...
                (2 * constants.F12 * ...
                superimposedParam.locStress(1, 1, jj) * ...
                superimposedParam.locStress(2, 1, jj)) + ...
                (constants.F66 * ...
                superimposedParam.locStress(3, 1, jj) ^ 2);
            
            if comparisonVar >= 1
                
                % Does laminate fail?
                fail{jj} = 'Yes';
               
                
            else
                
                % Does laminate fail?
                fail{jj} = 'No';
                
                % What is the mode of faliure?
                mode{jj} = 'N/A';
                
            end
        end
    end

end
