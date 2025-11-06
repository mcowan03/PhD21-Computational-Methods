%% Stationary Distribution for the Markov Process %%

% P' * x = x (where x is the invariant distribution)

function stationary = stationary_distribution(P)

[eigvecs, eigvals] = eig(P');

real_eigvals = real(diag(eigvals));


% eigenvector corresponding to largest real eigenvalue (eigvalue = 1)
%

[~, max_index] = max(real_eigvals);
stationary = eigvecs(:, max_index);

%ensure real
assert(all(imag(stationary) < 1e-10), 'Imaginary stationary distribution');
    stationary = real(stationary);
    
    % Normalize to sum to 1
    stationary = stationary / sum(stationary);
end