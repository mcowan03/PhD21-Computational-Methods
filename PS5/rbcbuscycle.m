%% AGGREGATE UNCERTAINTY %%

% TFP follows AR(1): log(A') = rho*log(A) +sigma_epsilon*epsilon
% where epsilon ~ N(0,1)

% Algorithm:
% 1. Parameterization
% 2. Simulate TFP path (for large T)
% 3. solution path for the model (n-th solution path) corresponding 
% to n-th state vector path (and so prices).
% 4. Solve backwards for expectation term in intertemporal equation, fill in counterfactuals
% 5. Given expectations, solve for optimal decision rules for all t.
% 6. Simulate forward
% 7. Check convergence, MSE between (n+1)-th and n-th price vector. If not, 
% repeat steps 4-7.

% We should load a steady state
% Then we look at the business cycle

% Params
sigma = 1;
rho = 0.95;
sigma_epsilon = 0.009;
alpha = 0.33;
beta = 0.99;
delta = 0.025;
chi = 1;
mu = 0.60;
eta = 7.8827;
muu = 0;

% Load steady state values
A_ss = 1;
[SS] = ssSolverNewA(A_ss);

% Aggregate TFP shock %
Nz = 7;
m=3;
% Tauchen %
[logzgrid, P] = tauchen(rho, sigma_epsilon, muu, Nz, m);
vGridA = exp(logzgrid)';
mTransA = P;


% Simulate path %
seed = 100;
rng(seed);
T = 10001;
BURNIN = 500; %initial burn-in period
requiredTime = T+BURNIN; %number of total periods
pInitialPoint = 1; 
vSimPath = fnSimulator(pInitialPoint,mTransA,BURNIN+T); % sim a path for TFP

vr      = SS.r*ones(requiredTime,1);
vw      = SS.w*ones(requiredTime,1);
vL      = SS.N*ones(requiredTime,1);
vC      = SS.C*ones(requiredTime,1);
vK      = SS.K*ones(requiredTime,1) + normrnd(0,0.0000001,requiredTime,1);
vY      = SS.Y*ones(requiredTime,1);
vI      = SS.K*delta*ones(requiredTime,1);

% separate paths to be updated iteratively
vCnew   = zeros(requiredTime,1);
vKnew   = SS.K*ones(requiredTime,1);

%=========================    
% preparation
%=========================    
% the updating weights
weightOld1   = 0.9900; % updating weight for capital stock 
weightOld2   = 0.9900; % updating weight for consumption

%=========================
% repeated transition method
%=========================
% iteration preparation    
iter    = 1;  % this is for interim reports
error2      = 10; % this is for terminal condition
errorTol = 1e-8;
maxIter = 10000;

% Vectorise shock-related paths
iA = vSimPath;
iFuture = [(2:requiredTime)';requiredTime];
iFutureFuture = [(3:requiredTime)';requiredTime;requiredTime];
futureShock = vSimPath(iFuture);
vA  = vGridA(iA); % Vector of A shocks

% prior calculation of time-series of the transition probabilities to the realized
% aggregate shocks on the simulated path
mTransRealized = zeros(size(vK));
for iTrans = 1:length(vK)
mTransRealized(iTrans,1) = mTransA(iA(iTrans),futureShock(iTrans));
end

% RTM loop
while error2 > errorTol && iter <= maxIter

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BACKWARDS SOLVE !
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % K next period - updated from new vK
    vKprime = [vK(2:end);vK(1)];
    
    % Empty object, tempV1 to carry cumsum expected values
    tempV1 = 0;

    for iAprime = 1:Nz
        Aprime = vGridA(iAprime);

        % Find a period where future shock realisation is same as iAprime
        % and capital stock is closest to vKpime from below/above
        % K value candidates for counterfactual, (K',A') states
        
        candidate = vK(find(vSimPath==iAprime)); % Recall vSimPath is index of A (7)
        candidateLocation = find(vSimPath==iAprime);
        candidate(candidateLocation>requiredTime-BURNIN) = []; % Last burnin periods cannot be candidates
        candidate(candidateLocation<BURNIN) = []; % Initial burnin periods cannot be candidates
        candidateLocation(candidateLocation>requiredTime-BURNIN) = []; % last burnin periods cannot be a candidate
        candidateLocation(candidateLocation<BURNIN) = [];  % initial burnin periods cannot be a candidate
        [candidate, index] = sort(candidate); % sort candidates in order
        candidateLocation = candidateLocation(index); % Save location

        % Now find period where capital stock is closest to vKprime from
        % below
        nLow = sum(repmat(candidate',length(vKprime),1)<vKprime,2); % Identifies all K candidates closest from below to Kprime
        nLow(nLow<=1) = 1; % the location cannot go below 1.
        nLow(nLow>=length(index)) = length(index)-1; % the location cannot go over the length(index)-1: note that it's the closest from BELOW.
        nHigh = nLow+1; %define the period where the capital stock is closest to vKprime from above
        weightLow = (candidate(nHigh) - vKprime)./(candidate(nHigh)-candidate(nLow)); %compute the weight on the lower side
        weightLow(weightLow<0) = 0; % optional restriction on the extrapolation
        weightLow(weightLow>1) = 1; % optional restriction on the extrapolation

        % Linearly interpolate beliefs (RHS of euler) given
        % counterfactual...
        % We want to calc wage, interest, labour participation

        vCLow = vC(candidateLocation(nLow));
        vCHigh = vC(candidateLocation(nHigh));
        vLLow = (((1-alpha).*Aprime.*vKprime.^alpha) ./ (eta.*vCLow.^sigma)).^(chi/(1+chi*alpha));
        vLHigh = (((1-alpha).*Aprime.*vKprime.^alpha) ./ (eta.*vCHigh.^sigma)).^(chi/(1+chi*alpha));
        vrLow = alpha.*Aprime.*(vKprime./vLLow).^(alpha-1) - delta;
        vrHigh = alpha.*Aprime.*(vKprime./vLHigh).^(alpha-1) - delta;

        psi2Low = -(mu/2)*((vKprime(candidateLocation(nLow ))./vKprime).^2-1);
        psi2High= -(mu/2)*((vKprime(candidateLocation(nHigh))./vKprime).^2-1);
        
        tempV1 = tempV1 + (futureShock ~= iAprime).* beta .*...
                mTransA(iA,iAprime).*...
                (weightLow.*(1./vC(candidateLocation(nLow ))).^(sigma).*(1+vrLow - psi2Low) ...
           + (1-weightLow).*(1./vC(candidateLocation(nHigh))).^(sigma).*(1+vrHigh - psi2High) );

    end

% Cumulatively update beliefs for actual realised future states
Lfuture    = ((1-alpha).*vGridA(futureShock).*vKprime.^alpha./(eta*vC(iFuture).^sigma)).^(chi/(1+alpha*chi));
rfuture    = alpha*vGridA(futureShock).*vKprime.^(alpha-1).*Lfuture.^(1-alpha) - delta;
psi2Future = -(mu/2)*((vKprime(iFuture)./vKprime).^2-1);
tempV1     = tempV1 + beta*...
                mTransRealized.*(1./vC(iFuture)).^(sigma).*(1+rfuture -  psi2Future); 
tempV1     = tempV1./(1+mu*(vKprime-vK)./vK);

% Update allocations
tempC = ((1./tempV1).^(1/sigma));
vL = ((1-alpha).*vA.*vK.^alpha./(eta*tempC.^sigma)).^(chi/(1+alpha*chi));
vw = (1-alpha).*vA.*vK.^(alpha).*vL.^(-alpha);
vr = alpha.*vA.*vK.^(alpha-1).*vL.^(1-alpha) - delta;
vI = (vr+delta).*vK + vw.*vL - tempC - (mu/2).*((vKprime-vK)./vK).^2.*vK;
vY = vA.*vK.^(alpha).*vL.^(1-alpha);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SIMULATE FORWARD
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        vKpast = [SS.K; vK(1:end-1)];
        
        vIpast = [SS.K*delta; vI(1:end-1)];
        vKnew = (1-delta)*vKpast + vIpast;
        vCnew = vGridA(vSimPath).*vKnew.^(alpha).*vL.^(1-alpha) - vI - (mu/2).*(vI./vKnew-delta).^2.*vKnew;
        
        % Calculate MSE
        error2 = mean(([vC - vCnew; vK - vKnew]).^2);
        errorK = vK - vKnew;

        % update
        vC      = weightOld1*vC + (1-weightOld1)*vCnew;
        vK      = weightOld2*vK + (1-weightOld2)*vKnew;


if (floor((iter-1)/50) == (iter-1)/50)
%=========================  
% Report
%========================= 
Phrase = ['Iteration is in progress: ',num2str(iter),'st iteration'];
disp(Phrase);
fprintf('Convergence criterion: \n');
fprintf('Error: %.18f \n', error2);
fprintf(' \n');
    
subplot(1,2,1);
plot(BURNIN:requiredTime-BURNIN,vK(BURNIN:requiredTime-BURNIN));hold on;
plot(BURNIN:requiredTime-BURNIN,vKnew(BURNIN:requiredTime-BURNIN),'-.');
xlim([BURNIN,requiredTime-BURNIN]);
hold off;
legend("Predicted K","Realized K","location","northeast");

subplot(1,2,2);
plot(BURNIN:requiredTime-BURNIN,vC(BURNIN:requiredTime-BURNIN));hold on;
plot(BURNIN:requiredTime-BURNIN,vCnew(BURNIN:requiredTime-BURNIN),'-.');
xlim([BURNIN,requiredTime-BURNIN]);
hold off;
legend("Predicted C","Realized C","location","northeast");

pause(0.2);

end

iter = iter+1;

end % end of the final loop
toc;

save '../rbcsolutions.mat';

%% Fitting LoM into linear specification %%

endoState = vK(BURNIN+1:end-1);
endoStatePrime = vK(BURNIN+2:end);
exoState = vGridA(vSimPath(BURNIN+1:end-1));
intercept = ones(size(endoState));

% Independent variable
x = [intercept, log(endoState), log(exoState), log(endoState).*log(exoState)];
% Dependent variable
y= log(endoStatePrime);

[coeff, bint1, r1, rint1, R1] = regress(y,x);
fprintf('======================== \n');
fprintf('Fitting the true LoM into the log-linear specification\n');
fprintf('======================== \n');
disp(['R-squared: ',num2str(R1(1))]);
fprintf(' \n');
fitlm(x,y,'Intercept',false)

%% Plot equilibrium capital path along with log-linear forecast (obtained above) %%

startingPoint = BURNIN+1;
startingEndo = endoState(1); % This is initial K
recovered = ones(1, requiredTime - BURNIN)*startingEndo;
for iTrans = 1:(requiredTime - BURNIN-1)
endoStateTemp = recovered(iTrans);
exoStateTemp = exoState(iTrans);

tempX = [1, log(endoStateTemp), log(exoStateTemp), log(endoStateTemp)*log(exoStateTemp)];
tempVal = coeff'*tempX';
recovered(iTrans+1) = exp(tempVal);
end

samplePeriod = 2500:10000;
figure;
plot(samplePeriod,endoState(samplePeriod),'Color','red','LineWidth',1.5);hold on;
plot(samplePeriod,recovered(samplePeriod),'Color','blue','LineStyle','--','LineWidth',1.5);
legend("True LoM","Linear LoM","location","best","FontSize",15)

%% Business cycle stats for raw time series %%
fprintf('\n');
fprintf('======================== \n');
fprintf('Business cycle statistics for the raw time series\n');
fprintf('======================== \n');
fprintf('mean log(output): %.4f \n', mean(log(vY)));
fprintf('st. dev. log(output): %.4f \n', std(log(vY)));
fprintf('skewness log(output): %.4f \n', skewness(log(vY)));
fprintf('------------------------ \n');
fprintf('mean log(investment): %.4f \n', mean(log(vI)));
fprintf('st. dev. log(investment): %.4f \n', std(log(vI)));
fprintf('skewness log(investment): %.4f \n', skewness(log(vI)));
fprintf('------------------------ \n');
fprintf('mean log(consumption): %.4f \n', mean(log(vC)));
fprintf('st. dev. log(consumption): %.4f \n', std(log(vC)));
fprintf('skewness log(consumption): %.4f \n', skewness(log(vC)));

%% Business cycle stats for HP filtered series %%
logvy = log(vY(BURNIN:requiredTime-BURNIN));
logvI = log(vI(BURNIN:requiredTime-BURNIN));
logvC = log(vC(BURNIN:requiredTime-BURNIN));
fprintf('\n');
fprintf('======================== \n');
fprintf('Business cycle statistics for the HP-filtered time series\n');
[~,vYhpfilter] = hpfilter(logvy); % Takes the cyclical component
[~,vIhpfilter] = hpfilter(logvI);
[~,vChpfilter] = hpfilter(logvC);
fprintf('mean log(output): %.4f \n', mean(log(vYhpfilter)));
fprintf('st. dev. log(output): %.4f \n', std(log(vYhpfilter)));
fprintf('skewness log(output): %.4f \n', skewness(log(vYhpfilter)));
fprintf('------------------------ \n');
fprintf('mean log(investment): %.4f \n', mean(log(vIhpfilter)));
fprintf('st. dev. log(investment): %.4f \n', std(log(vIhpfilter)));
fprintf('skewness log(investment): %.4f \n', skewness(log(vIhpfilter)));
fprintf('------------------------ \n');
fprintf('mean log(consumption): %.4f \n', mean(log(vChpfilter)));
fprintf('st. dev. log(consumption): %.4f \n', std(log(vChpfilter)));
fprintf('skewness log(consumption): %.4f \n', skewness(log(vChpfilter)));
