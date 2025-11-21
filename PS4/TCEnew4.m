%Aiyagari Perfect Foresight Transition (PFA)
% =========================================================================
% USING MORE INTUITIVE HISTOGRAM TRANSITION METHOD
% Solves for the transition path of aggregate capital (K) following an
% unannounced, permanent TFP shock from A=1.0 to A=1.1.
%
% ASSUMPTION: The following functions/variables are defined/available:
% - solve_aiyagari_srce(A_val): Solves the steady state for TFP=A_val.
% - tauchen(rho, sigma, mu, Nz, m): Generates logzgrid and P.
% - fnOptFAST(beta, budget, a_grid, expValue, Na, minWealth): Household optimization.
% =========================================================================
% THE FIX INVOLVED PROPERLY CALCULATING LABOUR SUPPLY !!! (rather than
% supplyL = 1)
% ======================================

% --- 0. PARAMETERS & GRIDS ---
alpha = 0.36;
delta = 0.08;
beta = 0.96;
a_lower = 0;
a_upper = 150;

% Grid Settings 
Na = 150;  % Coarse grid size for policy/value function (Backward Step)
Nz = 7;    % Productivity states
Na2 = 150; % Fine grid size for distribution (Forward Step)

% Non-linear grid setup 
x = linspace(0, 0.5, Na);
x2 = linspace(0, 0.5, Na2);
y = x.^5 / max(x.^5);
y2 = x2.^5 / max(x2.^5);
a_grid = a_lower + (a_upper - a_lower) * y';
a_grid2 = a_lower + (a_upper - a_lower) * y2';

% Productivity process 
rho = 0.9;
sigma = 0.2;
m = 3.0; 
mu = 0;
[logzgrid, P] = tauchen(rho, sigma, mu, Nz, m);
zgrid = exp(logzgrid)';
% Find stationary distribution, s.t. pi*P=pi
[eigvecs, eigvals] = eig(P'); % Returns eigenvalues and corresponding eigenvectors
[~, idx] = min(abs(diag(eigvals) - 1)); % find the eigenvalue equal (or close) to 1 (it = 1)
pi = eigvecs(:, idx); % corresponding eigenvector
pi = pi / sum(pi);                % normalize to sum to 1

supplyL = zgrid'*pi; %

% --- 1. SOLVE & LOAD INITIAL AND FINAL STEADY STATES ---
% We want:
% Distribution for original SRCE (A=1)
% New SRCE value function, V_T+1 (A=1.1)
T = 100;        % Transition horizon (periods)
error_path = 10;
tol_pfa = 1e-4;
gamma = 0.15;    % Damping parameter for updating the K path 
A_old = 1.0;
A_new = 1.1;

% Load Initial State (A=1.0) - check format of what we are loading.
fprintf('Solving for initial steady state (A=%.1f)...\n', A_old);
results_old = solve_aiyagari_srce2(A_old); % Unpack results.
K_old = results_old.K;
% Initial distribution (pi_old_flat) MUST be on the fine grid (Na2*Nz)
pi_old_flat_temp = results_old.currentDist0;
pi_old_flat_temp(pi_old_flat_temp < 0) = 0;
pi_old_flat = pi_old_flat_temp; % This should've already been done, but just in case
fprintf('Initial State solved: K_old = %.4f\n', K_old);

% Load Final State (A=1.1)
fprintf('Solving for final steady state (A=%.1f)...\n', A_new);
results_new = solve_aiyagari_srce2(A_new);
K_new = results_new.K;
mVF_new = results_new.mVF; % Terminal value function V_T (on Na x Nz grid)
fprintf('Final State solved: K_new = %.4f\n', K_new);

% --- 2. INDEXING AND PATH GUESS ---
numIterPFA = 1;

% Guess Initial K Path
K_path_old = zeros(T,1);
K_path_old(1) = K_old;
K_path_old(T) = K_new;
% Linear Guess 
K_path_old(2:T-1) = linspace(K_path_old(1), K_path_old(T), T-2)'; 

% --- 3. THE PERFECT FORESIGHT ALGORITHM LOOP ---
fprintf('\nStarting Perfect Foresight Algorithm (PFA) loop...\n');
while error_path > tol_pfa
    % 3a. Calculate Price Paths (r_t, w_t) based on the guessed K_t path
    r_path = alpha * A_new * (K_path_old / supplyL).^(alpha-1) - delta;
    w_path = (1 - alpha) * A_new * (K_path_old / supplyL).^alpha;
    
    % --- I. Household Optimization (Backward Induction) ---
    V_path = zeros(Na, Nz, T);
    aPrimePol_path = zeros(Na, Nz, T);
    
    V_path(:, :, T) = mVF_new; % Terminal condition V_T
    
    for t = T-1:-1:1
        r = r_path(t);
        w = w_path(t);
        V_future = V_path(:, :, t+1);
        
        for iz = 1:Nz
            z = zgrid(iz);
            % Expected future value: E[V(a', z')]
            expValue_mat = V_future * P';
            expValue = expValue_mat(:, iz); 
            minWealth = a_lower;
            
            for ia = 1:Na
                a = a_grid(ia);
                budget = w * z + (1 + r) * a;
                
                % Optimal Saving (using fnOptFAST)
                aprime = fnOptFAST(beta, budget, a_grid, expValue, Na, minWealth);
                aprime = max(aprime, a_lower); % Enforce borrowing constraint
                
                % Interpolate V(a') 
                ia_low = max(1, min(Na-1, sum(a_grid < aprime)));
                ia_high = ia_low + 1;
                denom = a_grid(ia_high) - a_grid(ia_low);
                
                if abs(denom) < 1e-10
                     weightLow = 1;
                else
                     weightLow = (a_grid(ia_high) - aprime) / denom;
                end
                weightLow = min(max(weightLow,0),1);
                weightHigh = 1 - weightLow;
                
                value = weightLow * expValue(ia_low) + weightHigh * expValue(ia_high);
                c = budget - aprime;
                
                if c <= 0
                    V_path(ia, iz, t) = -Inf; 
                    aPrimePol_path(ia, iz, t) = a_lower;
                else
                    V_path(ia, iz, t) = log(c) + beta * value;
                    aPrimePol_path(ia, iz, t) = aprime;
                end
            end 
        end
    end
       
    % --- II. Market Clearing (Forward Simulation) ---
    
    pi_path = zeros(Na2, Nz, T);
    K_path_new = zeros(T, 1);
    pi_path(:, :, 1) = pi_old_flat; % Start from the initial distribution
    
    for t = 1:T-1
        
        % Get current Na2*Nz distribution for current t:
        currentDist_t = pi_path(:, :, t);
        
        % Calculate realised capital at K_t
        marginalDista_t = sum(currentDist_t, 2);
        K_path_new(t) = sum(a_grid2(:) .* marginalDista_t(:));
        
        % 1. Interpolate Policy a'_t to finer grid (Na2)
        mAPrimePol_t_Na2 = zeros(Na2, Nz);
        for iz = 1:Nz
             mAPrimePol_t_Na2(:, iz) = interp1(a_grid, aPrimePol_path(:, iz, t), a_grid2, "linear", "extrap");
        end

        % Now calculate distribution at t+1 (we have current distribution
        % and policy function

        newDist_t_plus_1 = zeros(Na2, Nz); % Initialise an empty bucket for t+1
        for iz = 1:Nz
            for ia = 1:Na2
            % 1. Get the mass from the starting bucket (ia, iz) at time t
            mass = currentDist_t(ia, iz);
            
            % 2. Get the saving policy and interpolation weights
            aprime = mAPrimePol_t_Na2(ia, iz);
            
            LB = max(1, min(Na2-1, sum(a_grid2 < aprime)));
            UB = LB + 1;
            denom = a_grid2(UB) - a_grid2(LB);
            
            if abs(denom) < 1e-10
                weightLB = 1;
            else
                weightLB = (a_grid2(UB) - aprime) / denom;
            end
            weightLB = min(max(weightLB,0),1);
            weightUB = 1 - weightLB;
            
            % 3. Distribute this mass to the new buckets at t+1
                for iznext = 1:Nz
                    prob_z_trans = P(iz, iznext);

            % Add mass to the (LB, iznext) bucket
                    newDist_t_plus_1(LB, iznext) = newDist_t_plus_1(LB, iznext) + prob_z_trans * mass * weightLB;
                
                % Add mass to the (UB, iznext) bucket
                    newDist_t_plus_1(UB, iznext) = newDist_t_plus_1(UB, iznext) + prob_z_trans * mass * weightUB;
                end
            end
        end

        % Store new updated distribution in t+1:
        pi_path(:, :, t+1) = newDist_t_plus_1;

    end
    
    % Final K_T calculation
    currentDist_T = pi_path(:, :, T); 
    marginalDista_T = sum(currentDist_T, 2);
    K_path_new(T) = sum(a_grid2(:) .* marginalDista_T(:));
    
    % 3d. Check Convergence and Update
    error_path = max(abs(K_path_new - K_path_old));
    
    % Update the guess using damping (convex combination)
    K_path_old = gamma * K_path_new + (1 - gamma) * K_path_old;

    
    fprintf('PFA Iteration %d, Max Error: %.4e\n', numIterPFA, error_path);
    numIterPFA = numIterPFA + 1;
end

K_path_newSS = K_path_old;

% --- 4. PLOT TRANSITION PATH ---
figure;
plot(1:T, K_path_old, 'b-', 'LineWidth', 2);
hold on;
plot([1, T], [K_old, K_old], 'r--', 'LineWidth', 1);
plot([1, T], [K_new, K_new], 'g--', 'LineWidth', 1);
title(sprintf('Aiyagari Perfect Foresight Transition (A=%.1f to A=%.1f)', A_old, A_new));
xlabel('Time Period');
ylabel('Aggregate Capital (K_t)');
legend('K_t Path (Converged)', sprintf('Old SS K*=%.4f', K_old), sprintf('New SS K*=%.4f', K_new), 'Location', 'southeast');
grid on;
box on;
set(gca, 'FontSize', 12);
fprintf('\nTransition solved! Final error: %.4e\n', error_path);

% CALCULATING GINI COEFFICIENT %
% We have pi_path(Na, Nz, T).

gini_path = zeros(T, 1);
for ttt = 1:T
    azCurrentDist = pi_path(:, :, ttt); % Extracts Na*Nz distribution at time t
    azMarginalDist = sum(azCurrentDist, 2); % Extracts marginal distribution of wealth
    wealthMass_t = a_grid2 .* azMarginalDist; % Wealth mass at each level of a
    totalWealth_t = sum(wealthMass_t); % Should be = K_path_old_t
    % A check for this:
    % if ttt ==1 && totalWealth_t == K_path_old(t)
    % Get cumulative shares (add 0 for (0,0) origin)
    cum_pop_t = [0; cumsum(azMarginalDist)];
    cum_wealth_t = [0; cumsum(wealthMass_t)];
    
    % Normalize cumulative wealth to get Lorenz curve y-axis
    lorenz_y_t = cum_wealth_t / totalWealth_t;
    
    % Calculate area under the Lorenz curve (trapezoid rule)
    area_under_lorenz_t = trapz(cum_pop_t, lorenz_y_t);
    
    % 6. Calculate Gini and store it
    gini_path(ttt) = 1 - 2 * area_under_lorenz_t;
end
plot(1:T, gini_path, 'b-', 'LineWidth', 2);