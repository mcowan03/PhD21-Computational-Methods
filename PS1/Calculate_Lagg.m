function L_agg_output = Calculate_Lagg(w_input, T_input)

% --- Parameters (Keep these inside or make them global/passed as input) ---
a = 1;
alpha = 0.3;
tau = 0.15;
z_mean = 1;
A = 1;
r = 0.04;
beta = 0.96;

eta = 4.608;
chi = 0.956;
b = 0.062;
sigma = 0.276;

z_values = linspace(0.01, 10, 1000);
pdf_z = lognpdf(z_values, -0.5 * sigma^2, sigma);
pdf_z = pdf_z / sum(pdf_z); 

% --- Use Inputs ---
w = w_input; % This is the wage we are plotting against
T = T_input; % This is the fixed transfer for the curve

% --- Micro-level Solution (NO LOOP NEEDED) ---

% Store results
c_worker_values = zeros(size(z_values));
options = optimoptions('fsolve', 'Display', 'off');

% Loop through z values (Worker Intensive Margin)
for i = 1:length(z_values)
    z = z_values(i);
    % Define function handle for the equation f(c) = 0
    f = @(c) (1 + beta)*c - z * w * (1 - tau) * ((z * w * (1 - tau)) / (eta * c)^chi) - a * (1 + r * (1 - tau)) - T;
    c0 = 1;
    [c_sol, ~, exitflag] = fsolve(f, c0, options);
    
    if exitflag > 0 && c_sol > 0
        c_worker_values(i) = c_sol;
    else
        c_worker_values(i) = NaN;
    end
end

% Solve for labour hours
worker_hours = zeros(size(z_values));
for i = 1:length(worker_hours)
    % Only calculate if c is valid
    if ~isnan(c_worker_values(i)) && c_worker_values(i) > 0
        worker_hours(i) = ((z_values(i) * w * (1 - tau)) / (eta * c_worker_values(i)))^chi;
    else
        worker_hours(i) = NaN;
    end
end

% NON-WORKER SOLUTION & Value Function
% ... (Copy the c_non_worker, worker_value_func, and non_worker_value_func calculations) ...
% NOTE: You must include the check for NaN/invalid solutions in the worker_value_func as you did.
c_non_worker = zeros(size(z_values)); %initialise

    for i = 1:length(c_non_worker)
    c_non_worker(i) = (1 / (1 + beta)) * (b + a * (1 + r * (1 - tau)) + T);
    end

%disp(c_non_worker)

a_future_non_worker = beta * c_non_worker;

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



% Finding z* cutoff
tolerance = 0.005;
z_star_indices = find(abs(non_worker_value_func - worker_value_func) < tolerance);
if ~isempty(z_star_indices)
    z_star_index = min(z_star_indices);
else
    % Handle case where cutoff is outside grid (e.g., all work or all non-work)
    % For safety, assume minimum index (index 1) if no z* found within tolerance
    z_star_index = 1; 
end

% AGGREGATION
eff_labour_supply = worker_hours;
for i = 1:z_star_index
    eff_labour_supply(i) = 0; % Set labor to zero for z < z*
end

% Compute aggregate labor supply
L_agg = sum(eff_labour_supply .* pdf_z);
L_agg_output = L_agg;

end