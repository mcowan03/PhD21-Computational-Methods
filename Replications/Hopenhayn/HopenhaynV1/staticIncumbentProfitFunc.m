function [mStaticPi, nPol, mY] = staticIncumbentProfitFunc(eta, p , zgrid, w, cf)
% Calculates static labor demand, output, and profit.
% pi(z;p) = max_n{p * z * n^eta - w * n} - wcf

    % 1. Labor Policy (n): Analytically derived optimal labor choice
    % n(z) = ( (eta * p * z) / w )^(1 / (1 - eta))
    % Use element-wise operations (.*, .^, ./)
    nPol = ((eta * p .* zgrid) / w) .^ (1 / (1 - eta));
    
    % 2. Firm Output (Y)
    % Y(z) = z * n(z)^eta
    mY = zgrid .* nPol .^ eta;
    
    % 3. Static Profit (Pi)
    % Pi(z) = p*Y(z) - w*n(z) - cf
    mStaticPi = p .* mY - w .* nPol - w * cf;
end


