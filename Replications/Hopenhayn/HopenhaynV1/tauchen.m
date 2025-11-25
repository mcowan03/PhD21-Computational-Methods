%% TAUCHEN %%
% AR(1) process: x_t = (1 - rho) * mu + rho * x_{t-1} + epsilon_t, 
function [logzgrid, P] = tauchen(Nz, rho_z, sigma_z, z_bar, nStd)

% Discretise AR(1) = x(t+1) = mu + rho * x(t) + sigma*sqrt(1-rho^2)*e(t),
% Returns a vector: [logzgrid, P], where P is a transition matrix
% 1. Calculate Unconditional Standard Deviation (psi)
    % This is used to define the boundaries of the grid.
    psi = sigma_z / sqrt(1 - rho_z^2); 
    
    % 2. Calculate Innovation Standard Deviation (sigma_epsilon)
    % This is the standard deviation of the residual (epsilon_t)
    sigma_epsilon = sigma_z; % Assuming sigma is the standard deviation of the residual (epsilon)
    % NOTE: If sigma is defined as the unconditional std dev (as in some texts),
    % you would use: sigma_epsilon = sigma * sqrt(1 - rho^2); 
    % We assume the standard Tauchen interpretation where sigma is the shock std dev.
    
    % 3. Determine Step Size and Grid
    wstep = 2 * nStd * psi / (Nz - 1); % Step size (w)
    z_min = z_bar - nStd * psi;      % Minimum state
    
    logzgrid = linspace(z_min, z_min + (Nz - 1) * wstep, Nz);

    P = zeros(Nz, Nz);
    
    % 4. Calculate Transition Matrix (P)
    for i = 1:Nz % Current state x_t = Z_vals(i)
        
        % The expected mean of the next state (conditional mean)
        cond_mean = (1 - rho_z) * z_bar + rho_z * logzgrid(i);
        
        % Transition to the lowest state j=1
        % P(i, 1) = Prob(x_{t+1} <= Z_vals(1) + w/2)
        P(i, 1) = normcdf(logzgrid(1) + wstep/2, cond_mean, sigma_epsilon);
        
        % Transition to intermediate states j=2 to n-1
        % P(i, j) = Prob(Z_vals(j) - w/2 < x_{t+1} <= Z_vals(j) + w/2)
        for j = 2:(Nz-1)
            upper_bound = logzgrid(j) + wstep/2;
            lower_bound = logzgrid(j) - wstep/2;
            
            P(i, j) = normcdf(upper_bound, cond_mean, sigma_epsilon) - ...
                      normcdf(lower_bound, cond_mean, sigma_epsilon);
        end
        
        % Transition to the highest state j=n
        % P(i, n) = Prob(x_{t+1} > Z_vals(n) - w/2)
        P(i, Nz) = 1 - normcdf(logzgrid(Nz) - wstep/2, cond_mean, sigma_epsilon);
        
        % Normalize (needed for numerical stability)
        P(i, :) = P(i, :) / sum(P(i, :));
    end
end