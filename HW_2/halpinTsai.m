function [estimatedParameter] = halpinTsai(P, Pm, Pf, cf)

syms estimatedParameter

eta = (Pf - Pm) / (Pf - (estimatedParameter * Pm));

propEquation = P == (Pm * (1 + (estimatedParameter * eta * cf))) ./ ...
    (1 - (eta * cf));

estimatedParameter = double(solve(propEquation, estimatedParameter));

end
