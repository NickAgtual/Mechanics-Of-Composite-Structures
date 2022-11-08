function [deformationAtMidplane, z] = midplaneDeformation(...
    loading, Qbar, t)

% Z location of each ply
% Find way to automate this based on thickness
z = [-.254 -.127 0 .127 .254] * 10 ^ -3;

%% ABD Matrix

% Initializing A, B, and D matrices
[A, B, D] = deal(zeros(3, 3));

for ii = 1:3
    for jj = 1:3
        for kk = 2:length(t)+1
            
            A(ii, jj) = A(ii, jj) + (Qbar(ii, jj, kk-1) * (z(kk) - z(kk-1)));
            
            B(ii, jj) = B(ii, jj) + (.5 * (Qbar(ii, jj, kk-1) * ...
                (z(kk)^2 - z(kk-1)^2)));
            
            D(ii, jj) = D(ii, jj) + ((1/3) * (Qbar(ii, jj, kk-1) * ...
                (z(kk)^3 - z(kk-1)^3)));
            
        end
    end
end


%% Deflection caused by Uniaxial Load

% Laminate stiffness matrix
ABD = [A B;B D];

% Deflection due to uniaxial loading
deformationAtMidplane = inv(ABD) * loading';

end
