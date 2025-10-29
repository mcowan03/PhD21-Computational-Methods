%comp static for w* by varying A
% Parameterise
a = 1;
alpha = 0.3;
tau = 0.15;
z_mean = 1;
r = 0.04;
beta = 0.96;


% Want to estimate this via MSM, but to check code runs, set for now.
eta = 4.608;
chi = 0.956;
b = 0.062;
sigma = 0.276;

% Grid for z
z_values = linspace(0.01, 10, 1000);

% Parameterise for lognormal z distribution
mu = 0;
%sigma = 0.5;
pdf_z = lognpdf(z_values, -0.5 * sigma^2, sigma);
pdf_z = pdf_z / sum(pdf_z); % To normalise the sum of pdf_z(i)

% --- 1. Define A Grid ---
A_min = 0.5; 
A_max = 1.2;
num_A_points = 20; 
A_grid = linspace(A_min, A_max, num_A_points);

% --- 2. Initialize Storage ---
w_star_vector = zeros(size(A_grid));

% --- 3. Outer Loop: Solve GE for each A ---
for k = 1:length(A_grid)
    % Update the TFP parameter A
    A = A_grid(k);
    
   % Initialise guesses on prices and transfer
w = 0.2952;
T = -0.4292;

toleranceGE = 1e-5;
maxIterationsGE = 100;
iterationGE = 0;
error_L = 1;
error_G = 1;

%% THIS IS THE GE LOOP: Solves the model for guesses of w and T. Then checks errors at the aggregate.
% If there are errors - we adjust prices and cycle through the while loop.
while (abs(error_L) > toleranceGE || abs(error_G) > toleranceGE) && (iterationGE < maxIterationsGE)

% Store results
c_worker_values = zeros(size(z_values));

% Options for fsolve (suppress output)
options = optimoptions('fsolve', 'Display', 'off');

% Loop through z values
    for i = 1:length(z_values)
    z = z_values(i);

    % Define function handle for the equation f(c) = 0
    f = @(c) (1 + beta)*c - z * w * (1 - tau) * ((z * w * (1 - tau)) / (eta * c)^chi) - a * (1 + r * (1 - tau)) - T;

    % Initial guess for c (must be > 0)
    c0 = 1;

    % Use fsolve to find the root
    [c_sol, fval, exitflag] = fsolve(f, c0, options);

    % Store result if solution is valid and positive
        if exitflag > 0 && c_sol > 0
        c_worker_values(i) = c_sol;
        else
        c_worker_values(i) = NaN; % Mark invalid
        end
    end

% Display results
%disp(c_worker_values);

% Can solve for labour hours and savings
a_future = beta * c_worker_values;

%disp(a_future);

% Initialise
worker_hours = zeros(size(z_values));

    for i = 1:length(worker_hours)
    worker_hours(i) = ((z_values(i) * w * (1 - tau)) / (eta * c_worker_values(i)))^chi;
    end

%disp(worker_hours)

%% NON WORKER SOLUTION %%

c_non_worker = zeros(size(z_values)); %initialise

    for i = 1:length(c_non_worker)
    c_non_worker(i) = (1 / (1 + beta)) * (b + a * (1 + r * (1 - tau)) + T);
    end

%disp(c_non_worker)

a_future_non_worker = beta * c_non_worker;

%disp(a_future_non_worker)


%% Extensive margin decision %%

% Want to find where W(a,z) = N(a,z)

% Worker value function
worker_value_func = zeros(size(z_values));
for i = 1:length(worker_value_func)
    % Check for NaN or invalid solution from fsolve
    if isnan(c_worker_values(i)) || c_worker_values(i) <= 0 || isnan(worker_hours(i))
        % If the intensive margin is infeasible, utility from working is -infinity.
        worker_value_func(i) = -1e20; % Set to a large negative number to ensure W < N
    else
        % Original valid calculation
        a_future(i) = beta * c_worker_values(i);
        worker_value_func(i) = log(c_worker_values(i)) - eta * (1 / (1 + (1 / chi))) * worker_hours(i)^(1 + 1/chi) + beta * log(a_future(i));
    end
end

%disp(worker_value_func)

% Non-worker value function

non_worker_value_func = zeros(size(z_values));

    for i = 1:length(non_worker_value_func)
    non_worker_value_func(i) = log(c_non_worker(i)) + beta * log(a_future_non_worker(i));
    end

%disp(non_worker_value_func)

%% Finding z* cutoff
tolerance = 0.005

z_star = find(abs(non_worker_value_func - worker_value_func) < tolerance);
    if numel(z_star)>1
        z_star = min(z_star)
    else
        z_star = z_star;
    end


disp(z_star)

display(z_values(z_star))


%% AGGREGATION %%

% Construct a new labour supply vector, with the rule that h(a,z) = n(a,z)
% if z>=z*, 0 otherwise. Recall, we need the index which is == z_star

eff_labour_supply = worker_hours

    for i = 1:z_star
        eff_labour_supply(i) = 0;
    end
%disp(eff_labour_supply)

%% Aggregates

% Compute aggregate labor supply and aggregate capital
L_agg = sum(eff_labour_supply .* pdf_z)

a_vector = zeros(size(z_values));
a_vector(1:1000) = 1;

K_agg = sum(a_vector .* pdf_z) % = 1...


% Check error/discrepancy between aggregate supply and representative firm
% demand

% Firm demand for labour

L_demand = K_agg * (((1 - alpha) * A) / w)^(1 / alpha)

% Government Budget Balance
% Revenue
gov_revenue = tau * w * sum(eff_labour_supply .* z .* pdf_z) + a * r * tau
gov_spending = b * sum(pdf_z(1:z_star)) + T



% Compute errors
error_L = L_demand - L_agg;
error_G = gov_revenue - gov_spending;

    if abs(error_L) < toleranceGE && abs(error_G) < toleranceGE
    disp(['GE Found in ', num2str(iterationGE), ' iterations!'])
    disp(['Equilibrium Wage ', num2str(w)])
    disp(['Equilibrium Transfer ', num2str(T)])
    disp(['Equilibrium Labour Supply ', num2str(L_agg)])
    disp(['Equilibrium Labour Demand ', num2str(L_demand)])
    disp(['Government Spending ', num2str(gov_spending)])
    disp(['Government Revenue ', num2str(gov_revenue)])
        break;
    end

% Update prices

%compute implied wage and transfer

w_implied = (1 - alpha) * A * (K_agg / L_agg)^alpha;
step = 0.5;
w_new = w + step * (w_implied - w); 
T_new = T + error_G;

w = w_new;
T = T_new;

iterationGE = iterationGE + 1'

end
    w_star_vector(k) = w; % Store the converged equilibrium wage
end

% --- 6. Plot Comparative Statics ---
figure;
plot(A_grid, w_star_vector, 'k-', 'LineWidth', 2);
xlabel('Total Factor Productivity (A)');
ylabel('Equilibrium Wage (w*)');
title('Comparative Statics: Wage vs. Total Factor Productivity');
grid on;