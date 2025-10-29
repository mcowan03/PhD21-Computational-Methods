%% MSM SCRIPT %%

%Empirical moments
empirical_moments = [
    0.33;
    0.06;
    0.25;
    0.70
];

% Initial guess: eta, chi, b, sigma
params = [4.608, 0.956, 0.062, 0.276];

% Bounds: [eta, chi, b, sigma]
lb = [0.0, 0.0, 0.0, 0.0];
ub = [10, 10, 10, 10];

% Weighting matrix
W = diag([1, 1, 1, 1]);

% Define objective function handle
objective = @(params) msm_objective(params, empirical_moments, W);

% Run optimization
options = optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', 10);
[estimated_params, min_distance] = fmincon(objective, params, [], [], [], [], lb, ub, [], options);

disp('Estimated parameters:');
disp(estimated_params);
disp(simulate_model(estimated_params))
