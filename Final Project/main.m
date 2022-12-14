function [faliure, toExport, plotStress, plotStrain, z] = main(varargin)
%% Inputs 

% Structure fields (For lamina & strength properties)
inputStruct.fields = {'E1', 'E2', 'G12', 'v12', 'X', 'Xprime', 'Y', ...
    'Yprime', 'S', 'Xe', 'XePrime', 'Ye', 'YePrime', 'Se'};

% Structure conversion (for lamina and strength propertries)
inputStruct.conversion = [ones(1, 3) * 10 ^ 6, 1, ones(1, 5) * 10 ^ 3, ...
    ones(1, 5)];

% Loading identifiers
inputStruct.loading = {'Nx', 'Ny', 'Nxy', 'Mx', 'My', 'Mxy'};

% Initializing loading vecotr
loading = zeros(1, 6);

% If the code is executed via GUI
if nargin == 1
    
    % Boolean var signifying wether GUI was used
    GUI = 1;
    
    % Populating moduli and strength structures
    for ii = 1:length(inputStruct.fields)
        
        if (1 <= ii) && (ii <= 4)
    
        % Lamina properties
        moduli.(inputStruct.fields{ii}) = ...
            varargin{1}.(inputStruct.fields{ii}) * ...
            inputStruct.conversion(ii);
        
        elseif (5 <= ii) && (ii <= length(inputStruct.fields))
            
        % Strength properties
        strength.(inputStruct.fields{ii}) = ...
            varargin{1}.(inputStruct.fields{ii}) * ...
            inputStruct.conversion(ii);
        
        end
    end
    
    % Populating hygrothermal structure with local values
    hygrothermal.alpha = [varargin{1}.alpha1 varargin{1}.alpha2 0]; % 1/F
    
    % Populating hygrothermal structure with temperatures
    hygrothermal.T0 = varargin{1}.T0; % deg F
    hygrothermal.Tf = varargin{1}.Tf; % deg F
    
    % Reading loading from GUI input
    for ii = 1:length(inputStruct.loading)
        
        % Populating loading vector
        loading(ii) = varargin{1}.(inputStruct.loading{ii});
        
    end
    
    % Lamina thickness
    t = varargin{1}.PlyThickness;
    
    % Stacking sequence
    ss = varargin{1}.ss;
    
% If no inputs to main function, the defaul values are as follows
elseif nargin == 0
    
    % Boolean var signifying wether GUI was used
    GUI = 0;
    
    % Lamina Properties
    moduli.E1 = 7.4 * 10 ^ 6;
    moduli.E2 = 1.4 * 10 ^ 6;
    moduli.G12 = .73 * 10 ^ 6;
    moduli.v12 = .28;
    
    % Strength Properties
    strength.X = 260 * 10 ^ 3; % psi
    strength.Xprime = -162 * 10 ^ 3; % psi
    strength.Y = 6.2 * 10 ^ 3; % psi
    strength.Yprime = -32 * 10 ^ 3; % psi
    strength.S = 13.1 * 10 ^ 3; % psi
    strength.Xe = .03154;
    strength.XePrime = -.02189;
    strength.Ye = .00443;
    strength.YePrime = -.0229;
    strength.Se = .019;
    
    % Global Hygrothermal Properties
    hygrothermal.alpha = [3.9 * 10 ^ -6, 16.7 * 10 ^ -6, 0]; % 1/F
    
    % Temperature conditions
    hygrothermal.T0 = 250;
    hygrothermal.Tf = 100;
    
    % Loading
    loading = [0 10 0 5 0 0]; % 1:3 = forces 4:6 = moments
    
    % Lamina thickness
    t = .009;
    ss = [0 0 45 -45 -45 45 0 0]; % Stackup sequence

else
    
    % Error dialog box if invalid number of input arguments
    errordlg('Invalid Number of Input Arguments',...
        'Input Error');
    
end

%% Main Code Body

% Creating modified stackup sequence accounting for 'middle/zero' ply
[ssMod] = modifiedLayup(ss);

% Create Lamina strct with Qbar
[Qbar, S, Tepsilon] = reducedTransformedStiffnessMat(moduli, ...
    ssMod);

% Add Sbar to lamina struct
[Sbar, Tsigma] = reducedTransformedComplianceMat(S, ssMod);

% Create laminate strcture with deformationAtMidplane & ABD matrix
[deformationAtMidplane, z, ABD] = midplaneDeformation(loading, Qbar, ...
    ssMod, t, ss);

% Add global stress to lamina struct
[globLaminaStress] = globalLaminaStress(deformationAtMidplane, ...
    Qbar, z, ssMod);

% Add all remaining stress & strain to lamina struct
[laminaStressStrain] = transformation(globLaminaStress, Sbar, Tepsilon, ...
    Tsigma);

% Calculating hygrothermal stresses
[hygrothermal] = hygrothermalEffetcs(ssMod, hygrothermal, z, Tepsilon, ...
    Tsigma, ABD, Qbar, ss);

% Superimposing stress and strain
[superimposedParam] = superposition(laminaStressStrain, hygrothermal);

    % Plotting global stress and strain
[plotStress, plotStrain] = stressStrainPlots(superimposedParam, z, GUI);


% Checking faliure criteria
[faliure] = faliureCriteria(strength, superimposedParam);

% Effective mechanical and thermal laminate properties
[effectiveLaminatePropsMT] = effectiveLaminateProps(ABD, t, ss, ...
    hygrothermal);

% Arrange data for exporting
[toExport] = arrangeData(superimposedParam, z, deformationAtMidplane, ...
    hygrothermal, loading, inputStruct, effectiveLaminatePropsMT, GUI);

end