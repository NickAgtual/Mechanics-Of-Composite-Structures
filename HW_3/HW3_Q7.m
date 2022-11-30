%% HW3_Q7

%% Given Properties and Conditions
% Hygrothermal conditions
hygrothermal.T1 = 180; % deg C
hygrothermal.T2 = 30; % deg C

% Mechanical properties
lamina.E1 = 140; % GPa
lamina.E2 = 10; % GPa
lamina.G12 = 6; % GPa
lamina.v12 = .34;
lamina.t = 1 * 10 ^ -3; % m

% Hygrothermal properties
hygrothermal.alpha1 = 0;
hygrothermal.beta1 = 0;
hygrothermal.alpha2 = 30 * 10 ^ -6; % 1 / deg C
hygrothermal.beta2 = .55;

%% Reduced Transfromed Stiffness Matrix


