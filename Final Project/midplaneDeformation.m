function [deformationAtMidplane, z, ABD] = midplaneDeformation(...
    loading, Qbar, ssMod, t, ss, QbarT)

% Z location of each ply
if rem(length(ss), 2) == 0
    
    % Creating position of each lamina based on even # ply layup
    % Ex: [-2 -1 0 1 2]
    idx = [(-1 * flip(1:length(ss)/2)), 0, 1:length(ss)/2];
    
else
    
    % Creating position of each lamina based on even # ply layup
    % Ex: [-2 -1 -.5 .5 1 2]
    idx = [(-1 * flip(1:floor(length(ss)/2))), -.5, .5, ...
        1:floor(length(ss)/2)];
    
end

% z-coordinate
z = idx * t

%% ABD Matrix

% Initializing A, B, and D matrices
[A, B, D] = deal(zeros(3, 3));

for ii = 1:3
    for jj = 1:3
        for kk = 1:length(ssMod)
            
            if kk < (length(ss)/2) + 1
                
                A(ii, jj) = A(ii, jj) + (Qbar(ii, jj, kk) * ...
                    (z(kk + 1) - z(kk)));
                
                B(ii, jj) = B(ii, jj) + (.5 * (Qbar(ii, jj, kk) * ...
                    (z(kk + 1)^2 - z(kk)^2)));
                
                D(ii, jj) = D(ii, jj) + ((1/3) * (Qbar(ii, jj, kk) * ...
                    (z(kk + 1)^3 - z(kk)^3)));
                
            elseif kk == (length(ss)/2) + 1
                
                continue
                
            elseif kk > (length(ss)/2) + 1
                
                A(ii, jj) = A(ii, jj) + (Qbar(ii, jj, kk) * ...
                    (z(kk) - z(kk-1)));
                
                B(ii, jj) = B(ii, jj) + (.5 * (Qbar(ii, jj, kk) * ...
                    (z(kk)^2 - z(kk-1)^2)));
                
                D(ii, jj) = D(ii, jj) + ((1/3) * (Qbar(ii, jj, kk) * ...
                    (z(kk)^3 - z(kk-1)^3)));
            end
        end
    end
end


%% Deflection caused by Uniaxial Load

% Laminate stiffness matrix
ABD = [A B;B D];

% Deflection due to uniaxial loading
deformationAtMidplane = inv(ABD) * loading';

end
