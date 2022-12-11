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
            faliure(1).fail = 'Yes';
            
            % What is the mode of faliure?
            faliure(1).mode = 'Fiber Faliure';
            
        % Checking for matrix faliure (Method 1)
        elseif (min(superimposedParam.locStrain(1, 1, :)) ...
                <= strength.XePrime) || ...
                (min(supercapacitorParam.locStrain(1, 2, :)) >= ...
                strength.YePrie)
            
            % Does laminate fail?
            faliure(1).fail = 'Yes';
            
            % What is the mode of faliure?
            faliure(1).mode = 'Matrix ';
            
        elseif abs(max(superimposedParam.locStrain(1, 3, :))) >= ...
                strength.Se
            
            % Does laminate fail?
            faliure(1).fail = 'Yes';
            
            % What is the mode of faliure?
            faliure(1).mode = 'Matrix ';
            
        end 
    end

    function [out] = tsaiHillCriterion(strength, superimposedParam, ...
            faliure)
    
    end

    function [out] = tsaiWuCriterion(strength, superimposedParam, ...
            faliure)
        
    end

end
