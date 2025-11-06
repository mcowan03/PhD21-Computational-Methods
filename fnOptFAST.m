function [aprimeopt] = fnOptFAST(beta,budget,a_grid,expValue,Na,minWealth)

% optimization option
options = optimset('Display','off','TolFun',1e-12,'TolX',1e-12,'MaxIter',20000,'MaxFunEvals',20000); % MaxFunEvals: 22400 (Default);
[aprimeopt] = fminbnd(@residCal,minWealth,budget,options);

%=========================    
%nested function
%=========================    
function [optVal] = residCal(Guess)
    
aprime = Guess;
a_lower = sum(a_grid<aprime);
a_lower(a_lower<1) = 1;
a_lower(a_lower>=Na) = Na-1;
a_upper = a_lower+1;

weightLow = (a_grid(a_upper) - aprime)/(a_grid(a_upper)-a_grid(a_lower));
weightLow(weightLow<0) = 0;
weightLow(weightLow>1) = 1;

value = weightLow*expValue(a_lower)+(1-weightLow)*expValue(a_upper);
c = budget - aprime;

%update
value = log(c)+beta*value;
optVal = -value;

end    
%=========================    
%=========================    

end