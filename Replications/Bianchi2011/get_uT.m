function mu = get_uT(cT, yN, omega, eta, sigma)
    % Marginal Utility of Tradables
    % C = [w cT^-eta + (1-w) yN^-eta]^(-1/eta)
    term = omega * cT.^(-eta) + (1-omega) * yN.^(-eta);
    C = term.^(-1/eta);
    dC_dcT = C.^(1+eta) .* omega .* cT.^(-eta-1);
    mu = C.^(-sigma) .* dC_dcT;
end