function [aprimeopt] = fnOptFAST(beta,budget,a_grid,expValue,Na,minWealth)

% optimization option
options = optimset('Display','off','TolFun',1e-12,'TolX',1e-12,'MaxIter',20000,'MaxFunEvals',20000); % MaxFunEvals: 22400 (Default);
[aprimeopt] = fminbnd(@residCal, minWealth, budget, options);

    function [optVal] = residCal(Guess)
        aprime = Guess;
        
        ia_low = max(1, min(Na-1, sum(a_grid < aprime)));
        ia_high = ia_low + 1;

        denom = a_grid(ia_high) - a_grid(ia_low);
        if denom == 0
            weightLow = 1;
        else
            weightLow = (a_grid(ia_high) - aprime) / denom;
        end
        weightLow = min(max(weightLow, 0), 1);

        value = weightLow * expValue(ia_low) + (1-weightLow) * expValue(ia_high);
        c = budget - aprime;
        if c <= 0
            optVal = 1e6; % Large penalty if consumption non-positive
            return
        end

        value = log(c) + beta * value;
        optVal = -value;
    end

end