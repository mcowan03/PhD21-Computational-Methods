function distance = msm_objective(params, empirical_moments, W)
    % Simulate model moments for given params

   
    model_moments = simulate_model(params);

    % Moment errors
    moment_errors = model_moments - empirical_moments;

    % Weighting matrix W (identity or diagonal to start)

    W = diag([1, 1, 1, 1]);
    distance = moment_errors' * W * moment_errors;
end
