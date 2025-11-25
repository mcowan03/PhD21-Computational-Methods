%% VFI FOR INCUMBENT %%

function [VF, mStaticPi, mY, nPol, continuePol, exitCutoff] = vfiIncumbent(P, p, zgrid, beta, eta, w, cf, Nz, tol, maxit);

% Initialise

VFold = zeros(1, Nz);
VF = zeros(1, Nz);

% Solve static firm problem

[mStaticPi, nPol, mY] = staticIncumbentProfitFunc(eta, p , zgrid, w, cf);
% [logzgrid, P] = tauchen(Nz, rho_z, sigma_z, z_bar, nStd)

for it = 1:maxit
    expVFold = VFold * P';

    VF = mStaticPi + beta * max(expVFold, 0);

    dist = max(abs(VFold - VF));

    if dist < tol
        break
    end

    VFold = VF; % Update value function for the next iteration
end

% Continuing policy

EVF = VF * P';
continuePol = (EVF>=0); % binary vector: 1 if continue, 0 if exit

% Find productivity level corresponding to exit/entry
idx = find(continuePol, 1, 'first');
exitCutoff = zgrid(idx);

end