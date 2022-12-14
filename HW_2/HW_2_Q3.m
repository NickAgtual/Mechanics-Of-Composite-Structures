clear; close all; clc
%% Constituent Properties
prop.E2f = 14.8; % GPa
prop.Em = 3.45; % GPa
prop.vm = .36;

%% Transverse Modulus

% Range of fiber volume ratio
cf = linspace(0, 1, 100); % 100 values in range 0 - 1
% Matrix fiber volume
cm = 1 - cf;

% E'm 
EmPrime = prop.Em / (1 - (prop.vm) ^ 2);

% Transverse modulues (E2)
E2 = (prop.E2f * EmPrime) ./ ((cf .* EmPrime) + (cm .* prop.E2f));

%% Plotting

% Creating new figure
figure(1)

% Plotting E2/Em vs cf 
plot(cf, E2 ./ prop.Em, 'DisplayName', 'Carbon/Epoxy')

% Adding grid
grid on
grid minor

% Plot Descriptors
title('\emph{Transverse Modulus vs. Fiber Volume Ratio}', 'fontsize', ...
    16, 'Interpreter', 'Latex')
xlabel('\emph{$c_f$}', 'fontsize', 14, 'Interpreter', 'Latex')
ylabel('\emph{${\frac{E_{2f}}{E_m}}$}',...
    'fontsize', 14, 'Interpreter', 'Latex')
legend('location', 'best', 'Interpreter', 'Latex')

