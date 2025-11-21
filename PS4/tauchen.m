%% Tauchen method of discretising AR(1) process %%
function [logzgrid, P] = tauchen(rho, sigma, mu, Nz, m)

% Discretise AR(1) = x(t+1) = mu + rho * x(t) + sigma*sqrt(1-rho^2)*e(t), where e(t)~N(0,1)
%
% Arguments are the necessary parameters
%
% Returns a vector: [logzgrid, P], where P is a transition matrix




% Create a vector of z values, N points evenly spaced (7), with ub and lb
% defined by m std. deviations from mean. IS THIS WRONG??? sigmaz = 
logzgrid = linspace(mu - m * sqrt(sigma^2 / (1 - rho^2)), mu + m * sqrt(sigma^2 / (1 - rho^2)), Nz);

%z = linspace(mu - m * sqrt(sigma^2), mu + m * sqrt(sigma^2), N); % the
%alternative which I think may be correct...

%interval width:
z_step = logzgrid(2) - logzgrid(1);

% initialise transition matrix

P = zeros(Nz,Nz);

% Fill transition matrix

    for j = 1:Nz
        for k = 1:Nz
            if k == 1
            P(j, k) = normcdf((logzgrid(1) - mu - rho * logzgrid(j) + z_step/2)/ sigma);
            elseif k == Nz
            P(j, k) = 1 - normcdf((logzgrid(Nz) - mu - rho * logzgrid(j) - z_step / 2) / sigma);
            else
            P(j, k) = normcdf((logzgrid(k) - mu - rho * logzgrid(j) + z_step / 2) / sigma) - normcdf((logzgrid(k) - mu - rho * logzgrid(j) - z_step / 2) / sigma);
            end
        end
    end
end