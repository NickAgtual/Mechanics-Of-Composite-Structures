function [outputs] = faliureCriteria(strength, superimposedParam)

criteria = {'maxStressCriterion', 'maxStrainCriterion', ...
    'tsaiHillCriterion', 'tsaiWuCriterion'};


for ii = 1:length(criteria)
    
    faliure(1).criterion = criteria{ii};
    
end

    function [faliure] = maxStressCriterion(strength, ...
            superimposedParam, faliure)
        
        % Checking for fiber faliure
        if (max(superimposedParam.locStress(1, 1, :)) >= strength.X) || ...
                (max(supercapacitorParam.locStress(1, 2, :)) >= strength.Y)
            
            % Does laminate fail?
            faliure(1).fail = 'Yes';
            
            % What is the mode of faliure?
            faliure(1).mode = 'Fiber Faliure';
            
        % Checking for matrix faliure (Method 1)
        elseif (min(superimposedParam.locStress(1, 1, :)) ...
                <= strength.Xprime) || ...
                (min(supercapacitorParam.locStress(1, 2, :)) >= ...
                strength.Yprie)
            
            % Does laminate fail?
            faliure(1).fail = 'Yes';
            
            % What is the mode of faliure?
            faliure(1).mode = 'Matrix ';
            
        elseif abs(max(superimposedParam.locStress(1, 3, :))) >= strength.S
            
            % Does laminate fail?
            faliure(1).fail = 'Yes';
            
            % What is the mode of faliure?
            faliure(1).mode = 'Matrix ';
            
        end 
    end

    function [faliure] = maxStrainCriterion(strength, ...
            superimposedParam, faliure)
        
        % Checking for fiber faliure
        if (max(superimposedParam.locStrain(1, 1, :)) >= ...
                strength.Xe) || ...
                (max(supercapacitorParam.locStrain(1, 2, :)) >= ...
                strength.Ye)
            
            % Does laminate fail?
            faliure(2).fail = 'Yes';
            
            % What is the mode of faliure?
            faliure(2).mode = 'Fiber Faliure';
            
        % Checking for matrix faliure (Method 1)
        elseif (min(superimposedParam.locStrain(1, 1, :)) ...
                <= strength.XePrime) || ...
                (min(supercapacitorParam.locStrain(1, 2, :)) >= ...
                strength.YePrie)
            
            % Does laminate fail?
            faliure(2).fail = 'Yes';
            
            % What is the mode of faliure?
            faliure(2).mode = 'Matrix ';
            
        elseif abs(max(superimposedParam.locStrain(1, 3, :))) >= ...
                strength.Se
            
            % Does laminate fail?
            faliure(2).fail = 'Yes';
            
            % What is the mode of faliure?
            faliure(2).mode = 'Matrix ';
            
        end 
    end

    function [out] = tsaiHillCriterion(strength, superimposedParam, ...
            faliure)
    
    end

    function [out] = tsaiWuCriterion(strength, superimposedParam, ...
            faliure)
        
        % Tsai-Wu criterion constants
        constants.F11 = -1 / (strength.Xprime * strength.X);
        constants.F1 = (1 / strength.X) + (1 / strength.Xprime);
        constants.F22 = -1 / (strength.Yprime * strength.Y);
        constants.F2 = (1 / strength.Y) + (1 / strength.Yprime);
        constants.F12 = 1 / (2 * strength.Xprime * strength.X);
        constants.F66 = 1 / strength.S;
        
        comparisonVar = ...
            (constants.F1 * max(superimposed.locStress(1, 1, :))) + ...
            (constants.F2 * max(superimposed.locStress(1, 2, :))) + ...
            (constants.F11 * max(superimposed.locStress(1, 1, :))^ 2) + ...
            (constants.F22 * max(superimposed.locStress(1, 2, :))^ 2) + ...
            (2 * constants.F12 * max(superimposed.locStress(1, 1, :)) * ...
            max(superimposed.locStress(1, 2, :))) + ...
            (constants.F66 * max(superimposed.locStress(1, 3, :)) ^ 2);
        
        if comparisonVar >= 1
            
            % Does laminate fail?
            faliure(4).fail = 'Yes';
            
            % What is the mode of faliure?
            faliure(4).mode = 'Matrix ';
            
        end
            

    end

end
