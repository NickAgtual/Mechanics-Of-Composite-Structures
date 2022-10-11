%% Global System Definition

globSys.e1 = [1; 0 ; 0];
globSys.e2 = [0; 1; 0];
globSys.e3 = [0; 0; 1];

%% Local System Definition

syms theta

locSys.e1 = [0; 0; 1];
locSys.e2 = [1 / sqrt(2); 1 / sqrt(2); 0];
locSys.e3 = [1 / sqrt(2); -1 / sqrt(2); 0];

%% Transformation Matrix

globSys.comb = [globSys.e1 globSys.e2 globSys.e3];
locSys.comb = [locSys.e1 locSys.e2 locSys.e3];

for ii = 1: length(globSys.e1)
    for jj = 1: length(locSys.e1)
        
        b(ii, jj) = locSys.comb(jj, ii);
          
    end
end

%% Stress Transformation

globSys.stress = [3 2 1; 2 3 3; 1 3 5];

locSys.stress = b * globSys.stress * transpose(b);

