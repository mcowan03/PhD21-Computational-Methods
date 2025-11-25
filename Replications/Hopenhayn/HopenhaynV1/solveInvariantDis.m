%% INVARIANT DISTRIBUTION %%
% Solves for the stationary distribution of firms

function statDistribution = solveInvariantDis(P, pi_entrants, M, continuePol, Nz)
continuePol = continuePol(:);     % converts to column vector
F = diag(continuePol) * P;        % diag(continuepol) creates a 100x100 diag matrix with continuepol values on diag
% F is now a 100x100 matrix, with only positive transition values where the
% firm actually continues. (first rows will be 0 valued everywhere).

I = eye(Nz);
A = I - F';
B = M * pi_entrants(:);

statDistribution = A \ B;
statDistribution = statDistribution';
end