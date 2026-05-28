function [mu, d_mu] = get_uT_deriv(cT, yN, omega, eta, sigma)
    % Returns MU and numerical derivative d(MU)/dcT for Newton solver
    mu = get_uT(cT, yN, omega, eta, sigma);
    
    % Finite difference approximation for derivative
    eps = 1e-6;
    mu_up = get_uT(cT+eps, yN, omega, eta, sigma);
    d_mu = (mu_up - mu) / eps;
end