function val = get_uT_deriv_val(cT, yN, omega, eta, sigma)
    % Wrapper that returns ONLY the MU value for fzero
    % Also protects against negative inputs just in case
    if cT <= 0
        val = 1e10; % Return huge MU for negative C to push solver back
    else
        val = get_uT(cT, yN, omega, eta, sigma);
    end
end