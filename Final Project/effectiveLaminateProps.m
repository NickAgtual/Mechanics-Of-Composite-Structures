function [effectiveLaminateProps] = ...
    effectiveLaminateProps(ABD, t, ss, hygrothermal)

% Laminate total thickness
h = t * length(ss); % in

% abd matrix
abd = inv(ABD);

% Extracting a & b matrices
a = abd(1:3, 1:3);
b = abd(4:6, 1:3);

% Effective laminate properties (Mechanical)
effectiveLaminateProps.Ex = 1 / (h * a(1, 1));
effectiveLaminateProps.Ey = 1 / (h * a(2, 2));
effectiveLaminateProps.Gxy = 1 / (h * a(3, 3));
effectiveLaminateProps.vxy = -a(1, 2) / a(1, 1);

% Effective laminate properties (Thermal)
effectiveLaminateProps.therm = (a * hygrothermal.Nprelim) + ...
    (b * hygrothermal.Mprelim);

end
