%% Bianchi (2011) %%


% Parameters
r = 0.04;
sigma = 2;
eta = 1/0.83 - 1;
kappaN = 0.32;
kappaT = 0.32;
omega = 0.31;
beta = 0.91;

%% EXOGENOUS SHOCKS %%
% Shock process for yt and yn
Ny = 4;
rho_y = 0.8;
sigma_y = 0.058;
y_bar = 0;
nStd = 3;
% load("C:\Users\mc2500\OneDrive - University of Cambridge\PhD21 Computational Methods\Replications\Bianchi2011\shock_process.mat")
%yNGrid = yN;
%yTGrid = yT;
%P = shock_trans;
[logygrid, P] = tauchen(Ny, rho_y, sigma_y, y_bar, nStd);
yTGrid = exp(logygrid);
yNGrid = exp(logygrid);

% Simulation
T = 5000;
BURNIN = 0.1*T;
pathlength = T + BURNIN;
maxiter = 5000;
tol = 1e-5;
damp = 0.95;

%% SIMULATING A PATH FOR EXOGENOUS SHOCKS
iniPoint = 1;

yTSimPath = fnSimulator(iniPoint,P,pathlength);
yNSimPath = fnSimulator(iniPoint,P,pathlength);

%% Initial guess for allocation path
% Guess steady state
b0 = -0.2;
cTss = 1 + b0*r;
cNss = 1;

vB = b0*ones(pathlength, 1) + normrnd(0,0.0000001,pathlength,1); % Start with some debt
vC_T = cTss*ones(pathlength, 1);    % Initial guess for Consumption T
vC_N = ones(pathlength, 1); % Potentially unused
vMu = zeros(pathlength, 1); % Initial Mu (lagrange multiplier on occasionally binding constraint)

% Separate paths to be updated iteratively
vBnew = zeros(pathlength, 1);
vC_Tnew = zeros(pathlength, 1);
vMunew = zeros(pathlength, 1);


%% RTM LOOP %%
tic;
pNumIter = 0;
error = 2;

% Pre calculate
% Combine shocks to one index (independent shocks ~ easy: Prob(S'|S) = Prob(yT'|yT) * Prob(yN'|yN))
N_states = Ny*Ny;
vStatePath = (yTSimPath - 1) * Ny + yNSimPath; % Creates an index for all (yT, yN) pairs.
P_State = kron(P, P);

iS = vStatePath;
iFuture = [(2:pathlength)'; pathlength];
futureShock = vStatePath(iFuture);
vYN = yNGrid(yNSimPath)';
vYT = yTGrid(yTSimPath)';

% Pre-calculate Realized Transition Probabilities (for RTM exact match)
mTransRealized = zeros(pathlength, 1);
for t = 1:pathlength
    mTransRealized(t) = P_State(iS(t), futureShock(t));
end

% Prior calculation of time series of the transition probabilities to
% realised aggregate shocks on sim path


%

while error>1e-4



    % Backwards solution
    vMargUtil = get_uT(vC_T, vYN, omega, eta, sigma);
    vBprime = [vB(2:end); vB(1)];

    tempExpMU = zeros(pathlength, 1);

    % loop over possible future exogenous states S' = [yN', yT'], using the
    % index.

    for iStatePrime = 1:N_states

        % Find candidates where exogenous state was iStatePrime
        candidateLoc = find(vStatePath == iStatePrime);
        % Filter candidates (Burn-in handling similar to Lee's code)
        candidateLoc(candidateLoc > pathlength - BURNIN) = [];
        candidateLoc(candidateLoc < BURNIN) = [];
        if length(candidateLoc) < 2; continue; end % Hopefully not..

        % Find the endogenous state (B)
        candidate = vB(candidateLoc);
        
        % 2. Sort Candidates (to find closest neighbors)
        [candidate, index] = sort(candidate);
        candidateLoc = candidateLoc(index);

        BLow = sum(vBprime > candidate', 2); 
        
        % Boundary checks (from snippet)
        BLow(BLow <= 1) = 1;
        BLow(BLow >= length(index)) = length(index) - 1;
        BHigh = BLow + 1;
        % Weights interp
        candHigh = candidate(BHigh);
        candLow  = candidate(BLow);
        % weightLow = (K_High - K_target) / (K_High - K_Low)
        denom = candHigh - candLow;
        denom(denom == 0) = 1e-10; % Avoid div by zero
        weightLow = (candHigh - vBprime) ./ denom;
        % Clamp weights (from snippet)
        weightLow(weightLow < 0) = 0;
        weightLow(weightLow > 1) = 1;
        
        % 5. Interpolate MARGINAL UTILITY
        MU_Low  = vMargUtil(candidateLoc(BLow));
        MU_High = vMargUtil(candidateLoc(BHigh));
        MU_Interp = weightLow .* MU_Low + (1 - weightLow) .* MU_High;

        % Accumulate exp
        prob_trans = P_State(iS, iStatePrime);
        mask = (futureShock ~= iStatePrime); % Mask out realized future
        
        tempExpMU = tempExpMU + mask .* prob_trans .* MU_Interp;
    end

    % Future realisations
    MU_Realized = vMargUtil(iFuture); 
    tempExpMU = tempExpMU + mTransRealized .* MU_Realized;
    
    % EULER RHS: beta * (1+r) * E_t[u'(c_{t+1})]
    rhs_euler = beta * (1+r) * tempExpMU;

  % --- STEP 2: FORWARD SIMULATION (Notes Page 3, Step 5) ---
    % "Update state path & consumption... consistent with Euler Eq"
    for t = 1:pathlength
        yt = vYT(t);
        yn = vYN(t);
        bt = vB(t);
        rhs = rhs_euler(t);
        
        % Bounds for search (Prevents "Complex Function" crash)
        c_min = 1e-5;
        c_max = yt + (1+r)*max(abs(vB)) + 10; 

        % A. Unconstrained Solution (u'(c) = RHS)
        resid_unc = @(c) get_uT_deriv_val(c, yn, omega, eta, sigma) - rhs;
        
        try
            if resid_unc(c_min) * resid_unc(c_max) > 0
                c_unc = yt; % Root not bracketed
            else
                c_unc = fzero(resid_unc, [c_min, c_max]);
            end
        catch
            c_unc = yt; 
        end
        
        b_prime_unc = (1+r)*bt + yt - c_unc;
        
        % B. Check Constraint
        pN_unc = ((1-omega)/omega) * (c_unc/yn)^(eta+1);
        limit_val = -(kappaN * pN_unc * yn + kappaT * yt);
        
        if b_prime_unc >= limit_val
            % Not Binding
            vC_Tnew(t) = c_unc;
            vBnew(t+1) = b_prime_unc;
            vMunew(t)  = 0;
        else
            % C. Binding Constraint
            % Solve: (1+r)b + yT - c = -kappaN * pN(c)*yN - kappaT * yT
            resid_bind = @(c) ((1+r)*bt + yt - c) + ...
                              (kappaN * ((1-omega)/omega)*(c/yn)^(eta+1)*yn + kappaT * yt);
            try
                % Search strictly between 0 and c_unc
                if resid_bind(1e-6) * resid_bind(c_unc) > 0
                     % Fallback using fminbnd if sign doesn't change
                     c_bind = fminbnd(@(c) abs(resid_bind(c)), 1e-6, c_unc);
                else
                    c_bind = fzero(resid_bind, [1e-6, c_unc]);
                end
            catch
                c_bind = c_unc; 
            end
            
            vC_Tnew(t) = c_bind;
            vBnew(t+1) = (1+r)*bt + yt - c_bind;
            
            mu_val = get_uT_deriv_val(c_bind, yn, omega, eta, sigma) - rhs;
            vMunew(t) = max(0, mu_val);
        end
    end
    
    % Fix Indexing shift (vBnew(t+1) was calculated, restore vBnew(1))
    vBnew(1) = vB(1); 
    vBnew = vBnew(1:pathlength);
    
    % --- STEP 3: UPDATE (Notes Page 3, Step 6) ---
    diff = max(abs(vB - vBnew));
    error = diff;
    
    % Dampening
    vB   = damp*vB   + (1-damp)*vBnew;
    vC_T = damp*vC_T + (1-damp)*vC_Tnew;
    vMu  = damp*vMu  + (1-damp)*vMunew;
    
    if (floor((pNumIter-1)/50) == (pNumIter-1)/50)
        fprintf('Iter: %d | Error: %.6f\n', pNumIter, error);
        subplot(1,2,1);
        plot(1:pathlength,vB(1:pathlength));hold on;
        plot(1:pathlength,vBnew(1:pathlength),'-.');
        xlim([1,pathlength]);
        hold off;
        legend("Predicted B","Realized B","location","northeast");

        subplot(1,2,2);
        plot(1:pathlength,vC_T(1:pathlength));hold on;
        plot(1:pathlength,vC_Tnew(1:pathlength),'-.');
        xlim([1,pathlength]);
        hold off;
        legend("Predicted C","Realized C","location","northeast");
        
        pause(0.2);
    end
    
    pNumIter = pNumIter + 1;
    if pNumIter > maxiter; break; end
end
% --- End of Loop ---
toc;
%% PLOTTING RESULTS
figure;
subplot(3,1,1); plot(vB(BURNIN:end)); title('Bond Holdings'); grid on;
subplot(3,1,2); plot(vC_T(BURNIN:end)); title('Tradable Consumption'); grid on;
subplot(3,1,3); plot(vMu(BURNIN:end)); title('Crisis Multiplier (\mu)'); grid on;

%% policy function plots
%% BACKING OUT POLICY FUNCTIONS
% The RTM solution is a "cloud" of points. We visualize it by slicing
% the cloud for specific shock realizations.

% 1. Setup Data
b_today = vB(BURNIN:end-1);
b_tomorrow = vB(BURNIN+1:end);
c_today = vC_T(BURNIN:end-1);
state_idx = vStatePath(BURNIN:end-1);

% 2. Select Representative Shocks to Plot
% We will plot: Lowest Shock, Median Shock, Highest Shock
% (Based on diagonal yT=yN for simplicity)
diag_indices = find(yTGrid == yNGrid); % Find symmetric states
low_shock_idx = (1-1)*Ny + 1;          % Lowest (yT=1, yN=1)
mid_shock_idx = (2-1)*Ny + 2;          % Mid    (yT=2, yN=2)
high_shock_idx = (4-1)*Ny + 4;       % High   (yT=4, yN=4)

figure('Name', 'Backing Out Policy Functions');

% --- PLOT BOND POLICY b'(b) ---
subplot(1,2,1); hold on; grid on;
title('Bond Policy Function b''(b, y)');
xlabel('Current Debt b_t'); ylabel('Next Debt b_{t+1}');

% 45-degree line
plot([min(b_today), max(b_today)], [min(b_today), max(b_today)], 'w--', 'LineWidth', 1);

% Plot policies for specific shocks
plot_policy(b_today, b_tomorrow, state_idx, low_shock_idx, 'r', 'Low Shock');
plot_policy(b_today, b_tomorrow, state_idx, mid_shock_idx, 'g', 'Mid Shock');
plot_policy(b_today, b_tomorrow, state_idx, high_shock_idx, 'b', 'High Shock');
%legend('45-degree', 'Low (Crisis)', 'Mid', 'High', 'Location', 'NorthWest');

% --- PLOT CONSUMPTION POLICY c(b) ---
subplot(1,2,2); hold on; grid on;
title('Consumption Policy Function c(b, y)');
xlabel('Current Debt b_t'); ylabel('Consumption c_t');

plot_policy(b_today, c_today, state_idx, low_shock_idx, 'r', 'Low Shock');
plot_policy(b_today, c_today, state_idx, mid_shock_idx, 'g', 'Mid Shock');
plot_policy(b_today, c_today, state_idx, high_shock_idx, 'b', 'High Shock');

% --- Ergodic Debt Distribution --- %
%% PLOT ERGODIC DISTRIBUTION OF DEBT
figure('Name', 'Ergodic Distribution of Debt');

% 1. Remove Burn-in period (crucial for accuracy)
b_final = vB(BURNIN:end);

% 2. Histogram
subplot(2,1,1);
histogram(b_final, 50, 'Normalization', 'probability', 'FaceColor', [0.2 0.4 0.6]);
title('Histogram of Debt Holdings (b)');
xlabel('Bond Holdings (b)');
ylabel('Probability');
grid on;
xline(mean(b_final), 'r--', 'LineWidth', 2, 'Label', 'Mean');
% 3. Kernel Density (Smoothed version)
subplot(2,1,2);
[f, xi] = ksdensity(b_final);
plot(xi, f, 'LineWidth', 2, 'Color', 'k');
fill(xi, f, [0.8 0.8 0.8], 'FaceAlpha', 0.5); % Add shading
title('Smoothed Density of Debt Holdings');
xlabel('Bond Holdings (b)');
ylabel('Density');
grid on;

%% STATISTICS %%

%% CALCULATE CRISIS FREQUENCY
% Remove the burn-in period to calculate the true ergodic probability
mu_final = vMu(BURNIN:end);

% Define a small threshold to filter out numerical solver noise 
% (e.g., a solver might return mu = 1e-12 instead of exactly 0)
crisis_threshold = 1e-5; 

% Identify periods where the constraint binds
crisis_indicator = (mu_final > crisis_threshold);

% Calculate the frequency
crisis_frequency = sum(crisis_indicator) / length(mu_final);

fprintf('\n==================================================\n');
fprintf('Model Diagnostics:\n');
fprintf('Crisis Frequency: %.2f%%\n', crisis_frequency * 100);
fprintf('==================================================\n');


%% Consumption %%
%% CALCULATE AVERAGE CONSUMPTION
% Remove the burn-in period to calculate true ergodic moments
cT_final = vC_T(BURNIN:end);
cN_final = vYN(BURNIN:end); % Market clearing: cN = yN

% 1. Average Tradable Consumption
mean_cT = mean(cT_final);

% 2. Average Non-Tradable Consumption 
mean_cN = mean(cN_final);

% 3. Average Aggregate Consumption (CES Basket)
% Formula: c = [omega * (cT)^(-eta) + (1-omega) * (cN)^(-eta)]^(-1/eta)
c_basket_final = (omega * cT_final.^(-eta) + (1-omega) * cN_final.^(-eta)).^(-1/eta);
mean_c_basket = mean(c_basket_final);

%% CALCULATE CONSUMPTION VOLATILITY
% Calculate Standard Deviation (Volatility)
std_cT = std(cT_final);
std_cN = std(cN_final);
std_c_basket = std(c_basket_final);

%% CALCULATE AVERAGE DEBT
% Remove the burn-in period to calculate the true ergodic mean
b_final = vB(BURNIN:end);

% Calculate Average Debt
mean_b = mean(b_final);

% Print the result
fprintf('\n==================================================\n');
fprintf('Debt Diagnostics:\n');
fprintf('Average Debt (b): %.4f\n', mean_b);
fprintf('==================================================\n');

% Print the results
fprintf('\n==================================================\n');
fprintf('Consumption Volatility (Standard Deviation):\n');
fprintf('Volatility of Tradable Consumption (cT): %.4f\n', std_cT);
fprintf('Volatility of Non-Tradable Consumption (cN): %.4f\n', std_cN);
fprintf('Volatility of Aggregate Consumption (c): %.4f\n', std_c_basket);
fprintf('==================================================\n');
fprintf('\n==================================================\n');
fprintf('Consumption Diagnostics:\n');
fprintf('Average Tradable Consumption (cT): %.4f\n', mean_cT);
fprintf('Average Non-Tradable Consumption (cN): %.4f\n', mean_cN);
fprintf('Average Aggregate Consumption (c): %.4f\n', mean_c_basket);
fprintf('==================================================\n');

save('..\Bianchi2011\prelimresults.mat')


%% Checking KS1998 Algo for Bianchi (2011) %%
%==========================================================================  
% Fitting the aggregate Law of Motion (LoM) into a linear specification 
%========================================================================== 

% 1. Extract variables (removing burn-in)
endoState = vB(BURNIN+1:end-1);
exoStateT = vYT(BURNIN+1:end-1);
exoStateN = vYN(BURNIN+1:end-1);
endoStatePrime = vB(BURNIN+2:end);

intercept = ones(size(endoState));

% 2. Independent variables (LINEAR specification, no logs)
% x = [Constant, b_t, yT_t, yN_t, b_t*yT_t, b_t*yN_t]
x = [intercept, ...
     endoState, ...
     exoStateT, ...
     exoStateN, ...
     endoState .* exoStateT, ...
     endoState .* exoStateN];

% Dependent variable
y = endoStatePrime;

% 3. Run Regression
[coeff, bint1, r1, rint1, stats] = regress(y, x);
R_squared = stats(1);

fprintf('\n==================================================\n');
fprintf('Fitting the true LoM into the linear specification\n');
fprintf('==================================================\n');
disp(['R-squared: ', num2str(R_squared)]);
fprintf('\n');

% Optional: MATLAB's built-in table display for coefficients
fitlm(x, y, 'Intercept', false)

% 4. Recover the implied dynamics (Simulation using the estimated LoM)
simLength = length(endoState);
recovered = zeros(simLength, 1);

% Start at the true initial point
recovered(1) = endoState(1);

for iTrans = 1:(simLength - 1)
    b_temp  = recovered(iTrans);
    yT_temp = exoStateT(iTrans);
    yN_temp = exoStateN(iTrans);
    
    tempX = [1, ...
             b_temp, ...
             yT_temp, ...
             yN_temp, ...
             b_temp * yT_temp, ...
             b_temp * yN_temp];
             
    % Forecast next period's debt
    recovered(iTrans+1) = coeff' * tempX';
end

% 5. Plot the comparison
% Make sure sample period is within bounds
sample_start = min(500, floor(simLength/2));
sample_end   = min(1000, simLength);
samplePeriod = sample_start:sample_end;

figure('Name', 'Krusell-Smith LoM Check', 'Units', 'inches', 'Position', [1, 1, 10, 5]);
plot(samplePeriod, endoState(samplePeriod), 'Color', 'red', 'LineWidth', 1.5); hold on;
plot(samplePeriod, recovered(samplePeriod), 'Color', 'blue', 'LineStyle', '--', 'LineWidth', 1.5);

title('Aggregate Law of Motion: True vs Linear Forecast', 'FontSize', 14);
xlabel('Time (t)', 'FontSize', 12);
ylabel('Aggregate Debt (B)', 'FontSize', 12);
legend("True LoM", "Linear LoM", "Location", "best", "FontSize", 12);
grid on; hold off;

% Save figure (Ensure the folder exists)
if ~exist('../figures', 'dir')
    mkdir('../figures');
end
location = '../figures/lom_bianchi.pdf';
saveas(gcf, location);