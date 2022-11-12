%% HW3 Q4
clear; clc; close all

%% Material Properties

% Longitudinal modulus
lamina.E1 = 42.7 * 10 ^ 6; % psi

% Transverse moudulus
lamina.E2 = .92 * 10 ^ 6; % psi

% In-plane shear modulus
lamina.G12 = .71 * 10 ^ 6; % psi

% Poisson's ratio
lamina.v12 = .23;

%% 1 Calculate the A, B, D, a, b, d Matrices for the Layup [45/90/30/0]
