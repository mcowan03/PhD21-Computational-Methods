%% PS2 %%
% Discrete time, t=0 - t=50
% HH consumes, saves and supplies labour
% Human capital accumulates

% HH Decides whether to work or not at the start of a period
  % If unemployed at t-1, there is search effort to find a job, given by
  % phi
  % There is then a probability, psi, of getting a job.


% Parameterise
zeta = 0.01;
psi = 0.8;
phi = 0.5;
eta = 2;
b = 0.2;
h_bar = 10;
r = 0.04;
beta = 0.96;
sigma_z = 0.05;
rho = 0.9;
gamma_h = 0.1;
gamma_z = 0.1;
gamma_0 = 0.2;
econst = 0.577215;

% Discretize  idiosyncratic process for z by Tauchen Method
% 2 grid points on one std. dev range -> 2x2 transition matrix

%log(z) - N(0,sigma_z)

Nz = 2; %grid points
m = 1; % +- 1 std. dev.
sigma_x = sqrt(sigma_z)/sqrt(1-rho^2); %where x = log(z)
x = [-m*sigma_x, m*sigma_x];

zGrid = exp(x); % for levels, z. (zGrid(1) = exp(-sigma_x))

% Now need 2x2 transition matrix.

kappa = rho/sqrt(1-rho^2);
Phi = @(x) normcdf(x);

% Fill 2x2 transition matrix, P

P = zeros(Nz,Nz);
P(1,1) = Phi(kappa);
P(1,2) = 1 - P(1,1); % Probability of transitioning from state 1 to state 2
P(2,2) = P(1,1); % Probability of staying in state 2
P(2,1) = 1 - P(1,1); % Probability of transitioning from state 2 to state 1


% Set grids

% Wealth grid

pNumGrida = 30;
vGridamin = 0;
vGridamax = 30;
xx = linspace(0,0.5,pNumGrida);
yy = xx.^5/max(x.^5);
aGrid = vGridamin+(vGridamax-vGridamin)*yy;

% Human Capital grid
hGridPoints = 51;
hGrid = linspace(0,50, hGridPoints);

% Time
T=50;
num_periods = T+1;

%% Storage %%
% Must store solution for each state (a,h,z) at each time (t)

% --- Value Function Grids ---
% Initialize with NaN (Not-a-Number) to easily spot bugs
%4th dimension is for time
V = nan(pNumGrida, hGridPoints, Nz, num_periods);  % Value function if employed last period
S = nan(pNumGrida, hGridPoints, Nz, num_periods);  % Value function if unemployed last period

W = nan(pNumGrida, hGridPoints, Nz, num_periods);  % Choice-specific value of working
N = nan(pNumGrida, hGridPoints, Nz, num_periods);  % Choice-specific value of not working

% --- Policy Function Grids ---
a_prime_work_policy = nan(pNumGrida, hGridPoints, Nz, num_periods);
c_work_policy = nan(pNumGrida, hGridPoints, Nz, num_periods);

a_prime_not_work_policy = nan(pNumGrida, hGridPoints, Nz, num_periods);
c_not_work_policy = nan(pNumGrida, hGridPoints, Nz, num_periods);

% Wage function w(h,z)
wage_func = @(h, z) gamma_h * h + gamma_z * (z.^2) + gamma_0;

% Human Capital Accumulation funciton
h_prime_func = @(h) min(h + 1, h_bar);


% Problem at T = 50 %

t_index = num_periods; % T+1=51

[A, H, Z] = ndgrid(aGrid, hGrid, zGrid);

% Problem for worker:

%Recall: c + a' = w(h,z) + (1+r)a and euler...
income_W = wage_func(H,Z) + (1+r)*A;

%Optimal c and a':
a_prime_w = (beta / (1 + beta)) * income_W;
c_w = (1 / (1 + beta)) * income_W;

%calculate W_t
W_T_grid = log(c_w) - eta + beta * log(a_prime_w);

% Calculate Not Working problem

income_N = b + (1 + r) * A;

%Optimal c and a':
a_prime_n = (beta / (1 + beta)) * income_N;
c_n       = (1 / (1 + beta)) * income_N;

N_T_grid = log(c_n) + beta * log(a_prime_n);

% Using Gumbel Log-sum formula, we can calculate V_t and S_t:
V_T_grid = zeta * (econst + log(exp(W_T_grid / zeta) + exp(N_T_grid / zeta)));
v_W_S = (psi * W_T_grid + (1 - psi) * N_T_grid - phi); % for probabilities, psi.
S_T_grid = zeta * log(exp(v_W_S / zeta) + exp(N_T_grid / zeta));


% Store solutions

% Store the values
W(:, :, :, t_index) = W_T_grid;
N(:, :, :, t_index) = N_T_grid;
V(:, :, :, t_index) = V_T_grid;
S(:, :, :, t_index) = S_T_grid;

% Store the policies
a_prime_work_policy(:, :, :, t_index) = a_prime_w;
c_work_policy(:, :, :, t_index) = c_w;
a_prime_not_work_policy(:, :, :, t_index) = a_prime_n;
c_not_work_policy(:, :, :, t_index) = c_n;

%% LOOP %%

% --- Add necessary functions for interpolation and optimization ---

% Function 1: Interpolate V or S for the next period
% This uses linear interpolation across the a, h, and z dimensions.
% Since z' is discrete, we interpolate across a and h, and select the z' value.
% NOTE: The last argument to interpn must be a full grid of the input dimension.
% Since we pass a_prime (1x1) and h_prime (1x1), we pass the zGrid element too.
V_interp_func = @(a_prime, h_prime, next_period_value_grid, current_z_index) ...
    interpn(aGrid, hGrid, zGrid, next_period_value_grid, ...
            a_prime, h_prime, zGrid(current_z_index), 'linear');

% Function 2: Objective function for Not Working (N) interpolation
% We need a separate function for S since the job search probability (psi)
% and cost (phi) are included in the S_t calculation, not the N_t choice.
S_interp_func = @(a_prime, h_prime, next_period_value_grid, current_z_index) ...
    interpn(aGrid, hGrid, zGrid, next_period_value_grid, ...
            a_prime, h_prime, zGrid(current_z_index), 'linear');


%% 6. The Main Backward Loop

disp('Starting backward induction...');

% Loop backwards from t = T-1 down to t = 0
for t = T-1 : -1 : 0
    
    t_index = t + 1;       % Current time index (e.g., t=0 -> index 1)
    next_t_index = t + 2;  % Next time index (e.g., t=1 -> index 2)

    % Extract the NEXT period value grids for faster access
    V_next = V(:, :, :, next_t_index);
    S_next = S(:, :, :, next_t_index);
    
    % Loop through all current state grid points
    for i_a = 1:pNumGrida
        for i_h = 1:hGridPoints
            for i_z = 1:Nz
                
                % Current state values
                a = aGrid(i_a);
                h = hGrid(i_h);
                z = zGrid(i_z);
                
                % Expected human capital (deterministic)
                h_prime_w = h_prime_func(h); % If working
                h_prime_n = h;               % If not working
                
                % --- Step 1: Solve "Work" Problem (W_t) ---
                
                income_W = wage_func(h, z) + (1 + r) * a;
                
                % Define the Objective Function for Working (W)
                % Objective_W(a') = log(c) - eta + beta * E[V_t+1(a', h_prime_w, z')]
                objective_w = @(x) log(income_W - x) - eta + beta * ( ...
                    P(i_z, 1) * V_interp_func(x, h_prime_w, V_next, 1) + ...
                    P(i_z, 2) * V_interp_func(x, h_prime_w, V_next, 2) ...
                );
                
                % Find the maximization bounds for a'
                a_prime_lower_bound = vGridamin;
                % Upper bound ensures c > 0 (by a small epsilon)
                a_prime_upper_bound = income_W - 1e-6; 
                
                % Use fminbnd to find the a' that *maximizes* the objective.
                % MATLAB's fminbnd finds the MINIMUM, so we minimize the negative objective.
                [a_prime_opt_w, W_value] = fminbnd(@(x) -objective_w(x), ...
                                                   a_prime_lower_bound, a_prime_upper_bound);
                
                % Store optimal choice and value
                W(i_a, i_h, i_z, t_index) = -W_value; 
                a_prime_work_policy(i_a, i_h, i_z, t_index) = a_prime_opt_w;
                c_work_policy(i_a, i_h, i_z, t_index) = income_W - a_prime_opt_w;

                
                % --- Step 2: Solve "Not Work" Problem (N_t) ---
                
                income_N = b + (1 + r) * a;
                
                % Define the Objective Function for Not Working (N)
                % Objective_N(a') = log(c) + beta * E[S_t+1(a', h_prime_n, z')]
                objective_n = @(x) log(income_N - x) + beta * ( ...
                    P(i_z, 1) * S_interp_func(x, h_prime_n, S_next, 1) + ...
                    P(i_z, 2) * S_interp_func(x, h_prime_n, S_next, 2) ...
                );
                
                % Find the maximization bounds for a'
                a_prime_lower_bound = vGridamin;
                a_prime_upper_bound = income_N - 1e-6;
                
                [a_prime_opt_n, N_value] = fminbnd(@(x) -objective_n(x), ...
                                                   a_prime_lower_bound, a_prime_upper_bound);
                
                % Store optimal choice and value
                N(i_a, i_h, i_z, t_index) = -N_value; 
                a_prime_not_work_policy(i_a, i_h, i_z, t_index) = a_prime_opt_n;
                c_not_work_policy(i_a, i_h, i_z, t_index) = income_N - a_prime_opt_n;
                
                % --- Step 3: Calculate Ex-Ante Values (V_t, S_t) ---
                
                % V_t (Employed last period)
                % V_t = zeta * (econst + log( exp(W_t/zeta) + exp(N_t/zeta) ))
                V_t_value = zeta * (econst + log(exp(W(i_a, i_h, i_z, t_index) / zeta) + ...
                                                 exp(N(i_a, i_h, i_z, t_index) / zeta)));
                V(i_a, i_h, i_z, t_index) = V_t_value;


                % S_t (Unemployed last period)
                % The value of working for an S-state person is modified:
                v_W_S = (psi * W(i_a, i_h, i_z, t_index) + (1 - psi) * N(i_a, i_h, i_z, t_index) - phi);
                
                S_t_value = zeta * (econst + log(exp(v_W_S / zeta) + ...
                                                 exp(N(i_a, i_h, i_z, t_index) / zeta)));
                S(i_a, i_h, i_z, t_index) = S_t_value;

            end % End i_z loop
        end % End i_h loop
    end % End i_a loop
    
    disp(['Completed period t = ', num2str(t)]);
end % End time loop (t)