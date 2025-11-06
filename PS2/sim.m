%% Simulation %%

% --- Simulation Parameters ---
N_households = 20000; % Number of households to simulate
T = 50;               % Life cycle length

% Initial Conditions (t=0)
a0 = 1.0;             % Initial assets (example)
h0 = 0.0;             % Initial human capital
z0_index = 1;         % Assume all start in zGrid(1) (low shock state)

% --- Storage for Simulation Results ---
% Dimensions: (Time, Household)
Sim_A = nan(T + 1, N_households);
Sim_H = nan(T + 1, N_households);
Sim_Z = nan(T + 1, N_households);
Sim_L = nan(T + 1, N_households); % Labor supply (1=Work, 0=Not Work)
Sim_C = nan(T + 1, N_households); % Consumption

% Initialize state at t=0
Sim_A(1, :) = a0;
Sim_H(1, :) = h0;
Sim_Z(1, :) = z0_index;
Sim_L(1, :) = 1; % Assume all start employed (L_t-1 = 1) - this affects V/S choice at t=0

% Set up Random Number Streams for Shocks
rng('default'); % For reproducibility
Rand_Shock = rand(T, N_households);   % For Gumbel taste shocks (for labor choice)
Rand_Z_Draw = rand(T, N_households);  % For new productivity shock (z')

% You need a separate interpolation function for the policies, 
% since they are stored in separate grids
policy_interp_func = @(a, h, policy_grid, current_z_index) ...
    interpn(aGrid, hGrid, zGrid, policy_grid, a, h, zGrid(current_z_index), 'linear');

% You also need a function to calculate the probability of working (P_W)
% Remember P_W = exp(V_W/zeta) / (exp(V_W/zeta) + exp(V_N/zeta))
calc_prob_work = @(V_work, V_not_work, zeta) ...
    exp(V_work / zeta) ./ (exp(V_work / zeta) + exp(V_not_work / zeta));

% --- Main Simulation Loop (t=0 to T-1) ---

for t = 1:T % t=1 is the second row, t=0 is the first row
    t_model = t - 1; % Model time (0 to 49)
    t_index = t;     % MATLAB index for current period t
    
    % Next period index
    t_index_next = t + 1; 

    % Loop over all households
    for n = 1:N_households
        
        % Current State Variables
        a = Sim_A(t_index, n);
        h = Sim_H(t_index, n);
        z_idx = Sim_Z(t_index, n);
        L_prev = Sim_L(t_index, n); % Labor status in t-1 (1=Employed, 0=Unemployed)
        
        % 1. Determine Labor Supply Decision (L_t)
        
        % A. Look up the choice-specific values (W_t and N_t)
        W_t = W(i_a, i_h, i_z, t_index); % NOTE: MUST INTERPOLATE if a and h are continuous!
        N_t = N(i_a, i_h, i_z, t_index); % If you didn't force a/h on the grid
        
        % For simplicity here, we assume the simulated (a,h) are on the grid. 
        % In a real application, you must interpolate W_t and N_t! 
        
        % Find the nearest grid indices for a and h
        [~, i_a] = min(abs(aGrid - a));
        [~, i_h] = min(abs(hGrid - h));
        
        % Look up the W_t and N_t values directly from the stored grids
        W_t = W(i_a, i_h, z_idx, t_index); 
        N_t = N(i_a, i_h, z_idx, t_index); 

        % B. Calculate the value of working/not working based on previous status
        if L_prev == 1 % Employed last period (V-state)
            P_Work = calc_prob_work(W_t, N_t, zeta);
        else % Unemployed last period (S-state)
            % Value of working is modified (psi*W + (1-psi)*N - phi)
            V_Work_S = psi * W_t + (1 - psi) * N_t - phi;
            P_Work = calc_prob_work(V_Work_S, N_t, zeta);
        end
        
        % C. Draw the random shock and make the labor choice
        if Rand_Shock(t_model + 1, n) <= P_Work
            L_t = 1; % Work
        else
            L_t = 0; % Not Work
        end
        Sim_L(t_index_next, n) = L_t;
        
        
        % 2. Look up Policy Functions (a', c)
        
        % Choose the correct policy grid based on L_t
        if L_t == 1 % Work
            a_prime_grid = a_prime_work_policy(:, :, :, t_index);
            c_policy_grid = c_work_policy(:, :, :, t_index);
            h_prime = h_prime_func(h); % Human capital accumulates
        else % Not Work
            a_prime_grid = a_prime_not_work_policy(:, :, :, t_index);
            c_policy_grid = c_not_work_policy(:, :, :, t_index);
            h_prime = h; % Human capital stagnant
        end
        
        % Interpolate the optimal policies for the current continuous (a,h) state
        Sim_A(t_index_next, n) = policy_interp_func(a, h, a_prime_grid, z_idx);
        Sim_C(t_index_next, n) = policy_interp_func(a, h, c_policy_grid, z_idx);
        
        % 3. Determine Next State Variables (h', z')
        
        % Human Capital (deterministic)
        Sim_H(t_index_next, n) = h_prime; 
        
        % Productivity Shock (stochastic)
        % Get the transition probabilities for the current z_idx (Pi is P in your code)
        P_1_to_2 = P(z_idx, 2); % P(z' = z2 | z = z_idx)
        
        if Rand_Z_Draw(t_model + 1, n) <= P_1_to_2
            Sim_Z(t_index_next, n) = 2; % Transition to z2 (index 2)
        else
            Sim_Z(t_index_next, n) = z_idx; % Stay in z1 (index 1) OR transition to z1
        end
        
    end % End of household loop
end % End of time loop

