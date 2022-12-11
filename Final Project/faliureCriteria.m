function [outputs] = faliureCriteria(strength, superimposedParam)

criteria = {'maxStressCriterion', 'maxStrainCriterion', ...
    'tsaiHillCriterion', 'tsaiWuCriterion'};


for ii = 1:length(criteria)
    
    faliure(ii).criterion = criteria{ii};
    
end

    function [faliure] = maxStressCriterion(strength, ...
            superimposedParam, faliure)
        
        for jj = 1:length(superimposedParam)
            
            % Checking for fiber faliure
            if (superimposedParam.locStress(1, 1, jj) >= strength.X) || ...
                    (superimposedParam.locStress(1, 2, jj) >= strength.Y)
                
                % Does laminate fail?
                faliure(1).fail(jj) = 'Yes';
                
                % What is the mode of faliure?
                faliure(1).mode(jj) = 'Fiber Faliure';
                
                % Checking for matrix faliure (Method 1)
            elseif (superimposedParam.locStress(1, 1, jj) ...
                    <= strength.Xprime) || ...
                    (superimposedParam.locStress(1, 2, jj) >= ...
                    strength.Yprie)
                
                % Does laminate fail?
                faliure(1).fail(jj) = 'Yes';
                
                % What is the mode of faliure?
                faliure(1).mode(jj) = 'Matrix ';
                
            elseif abs(superimposedParam.locStress(1, 3, jj)) >= strength.S
                
                % Does laminate fail?
                faliure(1).fail(jj) = 'Yes';
                
                % What is the mode of faliure?
                faliure(1).mode(jj) = 'Matrix ';
                
            else
                
                % Does laminate fail?
                faliure(1).fail(jj) = 'No';
                
                % What is the mode of faliure?
                faliure(1).mode(jj) = 'N/A ';
            end
            
        end
    end

    function [faliure] = maxStrainCriterion(strength, ...
            superimposedParam, faliure)
        
        for jj = 1:length(superimposedParam)
            
            % Checking for fiber faliure
            if (superimposedParam.locStrain(1, 1, jj) >= ...
                    strength.Xe) || ...
                    (superimposedParam.locStrain(1, 2, jj) >= ...
                    strength.Ye)
                
                % Does laminate fail?
                faliure(2).fail(jj) = 'Yes';
                
                % What is the mode of faliure?
                faliure(2).mode(jj) = 'Fiber Faliure';
                
                % Checking for matrix faliure (Method 1)
            elseif (superimposedParam.locStrain(1, 1, jj) ...
                    <= strength.XePrime) || ...
                    (superimposedParam.locStrain(1, 2, jj) >= ...
                    strength.YePrie)
                
                % Does laminate fail?
                faliure(2).fail(jj) = 'Yes';
                
                % What is the mode of faliure?
                faliure(2).mode(jj) = 'Matrix ';
                
            elseif abs(superimposedParam.locStrain(1, 3, jj)) >= ...
                    strength.Se
                
                % Does laminate fail?
                faliure(2).fail(jj) = 'Yes';
                
                % What is the mode of faliure?
                faliure(2).mode(jj) = 'Matrix ';
                
            else
                
                % Does laminate fail?
                faliure(2).fail(jj) = 'No';
                
                % What is the mode of faliure?
                faliure(2).mode(jj) = 'N/A ';
                
            end
        end
    end

    function [out] = tsaiHillCriterion(strength, superimposedParam, ...
            faliure)
        
        for jj = 1:length(superimposedParam)
            
            comparisonVar = ...
                ((superimposedParam.locStress(1, 1, jj) / ...
                strength.X) ^ 2) + ...
                ((superimposedParam.locStress(1, 2, jj) / ...
                strength.Y) ^ 2) - ...
                ((superimposedParam.locStress(1, 1, jj) / strength.X) * ...
                (superimposedParam.locStress(1, 2, jj) / strength.X) + ...
                ((superimposedParam.locStress(1, 3, jj) / ...
                strength.S) ^ 2);
            
            if comprisonVar >= 1
                
                % Does laminate fail?
                faliure(3).fail(jj) = 'Yes';
                
                % What is the mode of faliure?
                faliure(3).mode(jj) = 'Matrix ';
                
            else
                
                % Does laminate fail?
                faliure(3).fail(jj) = 'No';
                
                % What is the mode of faliure?
                faliure(3).mode(jj) = 'N/A ';
                
            end
            
            
        end
        
    end

    function [out] = tsaiWuCriterion(strength, superimposedParam, ...
            faliure)
        
        for jj = 1:length(superimposedParam)
            % Tsai-Wu criterion constants
            constants.F11 = -1 / (strength.Xprime * strength.X);
            constants.F1 = (1 / strength.X) + (1 / strength.Xprime);
            constants.F22 = -1 / (strength.Yprime * strength.Y);
            constants.F2 = (1 / strength.Y) + (1 / strength.Yprime);
            constants.F12 = 1 / (2 * strength.Xprime * strength.X);
            constants.F66 = 1 / strength.S;
            
            comparisonVar = ...
                (constants.F1 * superimposed.locStress(1, 1, jj)) + ...
                (constants.F2 * superimposed.locStress(1, 2, jj)) + ...
                (constants.F11 * superimposed.locStress(1, 1, jj)^ 2) + ...
                (constants.F22 * superimposed.locStress(1, 2, jj)^ 2) + ...
                (2 * constants.F12 * superimposed.locStress(1, 1, jj) * ...
                superimposed.locStress(1, 2, jj)) + ...
                (constants.F66 * superimposed.locStress(1, 3, jj) ^ 2);
            
            if comparisonVar >= 1
                
                % Does laminate fail?
                faliure(4).fail(jj) = 'Yes';
                
                % What is the mode of faliure?
                faliure(4).mode(jj) = 'Matrix ';
                
            else
                
                % Does laminate fail?
                faliure(4).fail(jj) = 'No';
                
                % What is the mode of faliure?
                faliure(4).mode(jj) = 'N/A ';
                
            end
        end
    end

end
