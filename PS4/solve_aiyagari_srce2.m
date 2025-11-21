%% Aiyagari Script for different A %%
% Stores results.

function results = solve_aiyagari_srce2(A_value)
%% AIYAGARI 1994 SOLVER %%
verbose             = true; %reports 
A = A_value;
% Parameterize

    beta = 0.96;
    delta = 0.08;
    alpha = 0.36;
    

% Income process (z)
    m = 3.0; % m std. dev. either side of mu
    mu = 0;
    rho = 0.9;
    sigma = 0.2;
    Nz = 7; % number of grid points

%% Wealth Grid
    a_lower = 0.0;
    a_upper = 150.0;
    Na = 150;
    Na2 = 150;

    % Finer points at low wealth states: (Maliar, Maliar and Valli, 2010)
        x = linspace(0,0.5,Na);
        x2 = linspace(0,0.5,Na2);
        y = x.^5/max(x.^5);
        y2 = x2.^5/max(x2.^5);
        a_grid = a_lower+(a_upper-a_lower)*y;
        a_grid2 = a_lower+(a_upper-a_lower)*y2; %For interpolation

% Productivity
[logzgrid, P] = tauchen(rho, sigma, mu, Nz, m);
zgrid = exp(logzgrid);
    
    % Find stationary distribution, s.t. pi*P=pi
    [eigvecs, eigvals] = eig(P'); % Returns eigenvalues and corresponding eigenvectors
    [~, idx] = min(abs(diag(eigvals) - 1)); % find the eigenvalue equal (or close) to 1 (it = 1)
    pi = eigvecs(:, idx); % corresponding eigenvector
    pi = pi / sum(pi);                % normalize to sum to 1


        % Alternative method - repeated transition until pi=pi_next -
        % should be the same as the eigenvector method - sanity check
        % P is your 7x7 transition matrix
        tolmanual = 1e-12;
        maxitertransition = 1e5;

        % start with uniform row distribution
        pi_old = ones(1, size(P,1)) / size(P,1);

        for it = 1:maxitertransition
            pi_new = pi_old * P;                 % evolve forward
            diff = norm(pi_new - pi_old, 1);     % L1 distance
        if diff < tolmanual
            fprintf('Converged in %d iters, diff = %.2e\n', it, diff);
            break;
        end
        pi_old = pi_new;
        end

        if it == maxitertransition
            warning('Max iterations reached without convergence (diff = %.2e).', diff);
        end
  
  stationary = pi; % From eigenvector method.

%% INITIALISATION %%
  % Set up storage matrices
    mVF = zeros(Na, Nz); 
    mVFnew = mVF;
    mCPol = zeros(Na, Nz);
    mAPrimePol = zeros(Na, Nz);
    %define the size of the distribution
    currentDist = ones(Na2,Nz);
    currentDist = currentDist ./ (Na2*Nz); %normalize
  % Initial Guess
  % Increase speed for A=1.1
  if A == 1.1
      K = 8;
  else
      K = 7;
  end
    weightOld = 0.9;
    supplyL = zgrid*stationary;

%% SOLVING THE MODEL %%
% GE Loop Settings
error2 = 10;
tol_ge = 1e-5;
numIterGE = 1;

% GE LOOP %
while error2 > tol_ge
    % guessing K
    r = alpha*A*(K/supplyL)^(alpha-1)-delta;
    w = (1-alpha)*A*(K/supplyL)^alpha;

% VFI %
error = 10;
numIterVFI = 1;

    while error > 1e-5
        for iz = 1:Nz
               z = zgrid(iz);
               % expected future value, given z
               expValue = mVF * P';
               expValue = expValue(:, iz); % "column for this z"

               % We use monotonicity
               minWealth = a_lower;
                    
                   for ia = 1:Na

                   a = a_grid(ia);
                   budget = w*z+(1+r)*a;

                   % Optimal Saving
                   aprime = fnOptFAST(beta,budget,a_grid,expValue,Na,minWealth);

                   % Borrowing constraint
                   if aprime < a_lower
                       aprime = a_lower;
                   end

                   % Linear interpolation for off-the-grid aprime
                   ia_low = max(1, min(Na-1, sum(a_grid < aprime)));
                   ia_high = ia_low + 1;
                   denom = a_grid(ia_high) - a_grid(ia_low);
                   if denom == 0
                       weightLow = 1;
                   else
                       weightLow = (a_grid(ia_high) - aprime) / denom;
                   end
                   weightLow = min(max(weightLow,0),1);

                   value = weightLow*expValue(ia_low)+(1-weightLow)*expValue(ia_high);

                   c = budget - aprime;

                   % Update
                   minWealth = aprime;
                   mVFnew(ia,iz) = log(c)+beta*value;
                   mCPol(ia,iz) = c;
                   mAPrimePol(ia,iz) = aprime;

                   end
        end

% Iteration
error = max(max(max(abs(mVFnew-mVF))));
mVF = mVFnew;
numIterVFI = numIterVFI + 1;

    end % Ends VFI loop

if verbose
    fprintf('VFI converged after %d iterations (final error: %.2e)\n', numIterVFI-1, error);
end


 %interpolation - wealth policy
    %=========================
    %define the size of interpolated policy function
    mAPrimePol2 = zeros(size(currentDist));
    %interpolate
    if Na == Na2
        mAPrimePol2 = mAPrimePol;
        else
        for iz = 1:Nz
        mAPrimePol2(:,iz) = interp1(a_grid,mAPrimePol(:,iz),a_grid2,"linear","extrap");
        end
    end

% Non Stochastic Simulation - Iterative Method % Similar to KFE?
% Going to find a large transition probability matrix %

errDist = 20;
numIterDist = 1;
pTolDist = 1e-7;
while errDist > pTolDist

    mNewDist = zeros(size(currentDist)); % Reset distribution, we are going to move mass into this new distribution

    for izz = 1:Nz
        for iaa = 1:Na2
            mass = currentDist(iaa, izz); % Get mass at (iaa, izz)

            aprime = mAPrimePol2(iaa,izz); % Gets new a'
            % Interpolation and weights
            LB = max(1, min(Na2-1, sum(a_grid2 < aprime))); %index
            UB = LB + 1;
            denom = a_grid2(UB) - a_grid2(LB);
            if abs(denom) < 1e-10 % Use robust check
                weightLB = 1;
            else
                weightLB = (a_grid2(UB) - aprime) / denom;
            end
            weightLB = min(max(weightLB,0),1);
            weightUB = 1 - weightLB;

            % Distribute mass to new states
            for iznext = 1:Nz
                prob_z_trans = P(izz, iznext); % Prob of z' given z

                mNewDist(LB, iznext) = mNewDist(LB, iznext) + prob_z_trans*mass*weightLB;
                mNewDist(UB, iznext) = mNewDist(UB, iznext) + prob_z_trans*mass*weightUB;
            end
        end
    end

    errDist = max(max(abs(mNewDist - currentDist)));
    currentDist = mNewDist; % Update for next iteration
    numIterDist = numIterDist + 1;
end


%endogenous aggregate allocation
marginalDista = sum(currentDist,2); %  Sums each element in a row, to get mass for a particular asset level. (gives column vector of weights at each asset level).
endoK = a_grid2*marginalDista;

%error and update
error2 = abs(endoK - K);
K = K.*weightOld+endoK.*(1-weightOld);

% report only spasmodically
if verbose == true && (floor((numIterGE-1)/50) == (numIterGE-1)/50) || error2<= tol_ge
%=========================  
% interim report
%=========================  

fprintf(' \n');
fprintf('market clearing results \n');
fprintf('max error: %.15f \n', error2);
fprintf('capital rent: %.15f \n', r);
fprintf('wage: %.15f \n', w);
fprintf('aggregate capital: %.15f \n', K);

% plot
close all;
figure;
plot(a_grid2,currentDist);
title("The wealth distributions for different labor endowments","fontsize",15)
hold off;

pause(0.01);
toc;

end

numIterGE = numIterGE + 1;

end


% 3. Output Results
    % Collect the necessary final, converged values
    results.K = K;                       % Aggregate Capital (K*)
    results.r = r;                       % Interest Rate (r*)
    results.w = w;                       % Wage (w*)
    results.mVF = mVF;                   % Final Value Function V*(a, z)
    results.mAPrimePol = mAPrimePol;     % Policy Function a'(a, z) (Na x Nz)
    
    currentDist0 = currentDist; % is this actually necessary...

    results.currentDist0 = currentDist0; % Flattened Stationary Distribution pi*
    results.mAPrimePol2 = mAPrimePol2;   % Policy Function on fine grid (Na2 x Nz)
    
end