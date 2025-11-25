%% HOPENHAYN 1992 %%


%% Parameterise %%
beta = 0.8;
eta = 2/3;
cf = 20;
ce = 40;
D = 100;
w = 1;
rho_z = 0.9;
sigma_z = 0.2;
Nz = 100;
z_bar = 1.3;
nStd = 4;

M = 1;
% Tolerances/Settings

tol = 1e-7;
maxit = 200;

%% Discretise TFP process
[logzgrid, P] = tauchen(Nz, rho_z, sigma_z, z_bar, nStd);
zgrid = exp(logzgrid);


%% Initial Productivity for Entrants
% Use stationary distribution from our AR(1) process... (tbd)
[eigvecs, eigvals] = eig(P'); % Returns eigenvalues and corresponding eigenvectors
[~, index] = min(abs(diag(eigvals) - 1)); % find the eigenvalue equal (or close) to 1 (it = 1)
pi_entrants = eigvecs(:, index); % corresponding eigenvector
pi_entrants = pi_entrants / sum(pi_entrants);                % normalize to sum to 1


%% Stationary Equilibrium %%
ptol = 1e-5;
pmin = 1;
pmax = 100;
iterp = 0;

while (pmax - pmin) > ptol && iterp < maxit
    iterp = iterp + 1;
    % Guess a price:
    p = (pmin + pmax)/2;

    % Calculate incumbent VF
    [VF, ~, ~, ~, ~, ~] = vfiIncumbent(P, p, zgrid, beta, eta, w, cf, Nz, tol, maxit);

    % Calculate entrant firm value 
    VFentrant = beta * (VF*pi_entrants);

    % Check free entry condition (FEC)
    diff = VFentrant - ce;
    if abs(diff) < ptol
        break; % Price converged
    elseif diff > 0
        % VFentrant > ce. Entry is too profitable -> Price needs to fall (pmax = price)
        pmax = p;
    else
        % VFentrant < ce. Entry is not profitable -> Price needs to rise (pmin = price)
        pmin = p;
    end
end

if iterp == maxit
    warning('Price bisection did not converge within %d iterations.', maxit);
end


pss = p;
fprintf('Bisection Converged in %d iterations.\n', iterp);
fprintf('Equilibrium Price (p*): %.4f\n', pss);

%% Final Equilibrium Calculations (using p*) %%
[VF, firm_profit, firm_output, pol_n, pol_continue, exit_cutoff] = vfiIncumbent(P, pss, zgrid, beta, eta, w, cf, Nz, tol, maxit);
% Stationary distribution evolution
Minitial = 1;
distrib_stationary_0 = solveInvariantDis(P, pi_entrants, Minitial, pol_continue, Nz);
mStar = D / (distrib_stationary_0 * firm_output');

% Rescale
distrib_stationary = mStar * distrib_stationary_0;
totalMass = sum(distrib_stationary);

% Normalise for percent
pdfStationary = distrib_stationary / totalMass;
cdfStationary = cumsum(pdfStationary);


% Solving for M using D(p) = Y(p,M)

%% 4. EMPLOYMENT DISTRIBUTIONS %%
% f. calculate employment distributions
% pol_n (labor demand) * distrib_stationary (mass of firms) gives total labor demand in each state.
distrib_emp = pol_n .* distrib_stationary;
 
% invariant employment distribution by percent
pdf_emp = distrib_emp / sum(distrib_emp);
cdf_emp = cumsum(pdf_emp);

%% 5. FINAL STATISTICS AND OUTPUT %%
total_employment = pol_n * distrib_stationary';
average_firm_size = total_employment / totalMass;
exit_rate = mStar / totalMass; % Entry rate equals exit rate in steady state (m*/total_mass)


fprintf('\n-----------------------------------------\n');
fprintf('Stationary Equilibrium Results\n');
fprintf('-----------------------------------------\n');
fprintf('Equilibrium Price (p*): %.4f\n', pss);
fprintf('Entry/Exit Rate: %.3f\n', exit_rate);
fprintf('Total Mass of Firms: %.2f\n', totalMass);
fprintf('Average Firm Size: %.2f\n', average_firm_size);
fprintf('Exit Cutoff (Productivity): %.2f\n', exit_cutoff); 
fprintf('Total Employment (Demand): %.2f\n', total_employment);


%% PLOTTING %%


figure;
    
    % 1. Plot V(z)
    plot(zgrid, VF, 'LineWidth', 2, 'Color', '#0072BD', 'DisplayName', 'Value Function'); % Blue line
    hold on;
    
    % 2. Plot exit threshold
    xline(exit_cutoff, 'r--', 'LineWidth', 1.5, 'Alpha', 0.7, ...
          'DisplayName', ['Exit Threshold z^* = ' num2str(round(exit_cutoff, 2))]);
      
    % 3. Zero plot
    % axhline equivalent
    yline(0, 'k--', 'LineWidth', 1, 'Alpha', 0.6, 'DisplayName', 'Zero Value'); % Black dashed line
    
    % Add Title, Labels, and Legend
    title('Incumbent Firm Value Function V(z)');
    xlabel('Productivity Level (z)');
    ylabel('Value (Present Discounted Profit)');
    legend('Location', 'best');
    grid on;
    
    hold off;

% --- 2. Stationary PDF (Firms vs. Employment) ---
    figure;
    
    % 1. Plot the Firm Distribution PDF 
    plot(zgrid, pdfStationary, 'LineWidth', 2, 'Color', '#0072BD', 'DisplayName', 'Share of Firms');
    hold on; 
    
    % 2. Plot the Employment Distribution PDF 
    plot(zgrid, pdf_emp, 'LineWidth', 2, 'Color', '#D95319', 'DisplayName', 'Share of Employment');
    
    title('Stationary Distributions: Firms vs. Employment (PDF)');
    xlabel('Productivity Level');
    ylabel('Probability Density');
    legend('Location', 'best');
    grid on;
    hold off; 
 
    %% --- PIE CHART CALCULATIONS  ---

employed_cutoffs = [20, 50, 100, 500];
num_bins = length(employed_cutoffs) + 1;
share_firms = zeros(num_bins, 1);
share_employment = zeros(num_bins, 1);
labels = {'<20','21-50','51-100','101-500','>500'};

% Prepare data for interpolation: pol_n must be sorted and unique
[sorted_n, I] = unique(pol_n);
cdf_stat_sorted = cdfStationary(I);
cdf_emp_sorted = cdf_emp(I);

% --- Calculate Shares by Interpolating CDFs ---

% Firm Shares
previous_cdf_value_firm = 0;
for i = 1:length(employed_cutoffs)
    cutoff = employed_cutoffs(i);
    % Interpolate the CDF of firms at the cutoff labor demand
    interpolate = interp1(sorted_n, cdf_stat_sorted, cutoff, 'linear', 'extrap');
    
    current_share = interpolate - previous_cdf_value_firm;
    share_firms(i) = max(0, current_share);
    previous_cdf_value_firm = previous_cdf_value_firm + current_share;
end
share_firms(num_bins) = max(0, 1 - previous_cdf_value_firm);
share_firms = share_firms / sum(share_firms); % Re-normalize

% Employment Shares
previous_cdf_value_emp = 0;
for i = 1:length(employed_cutoffs)
    cutoff = employed_cutoffs(i);
    % Interpolate the CDF of employment at the cutoff labor demand
    interpolate = interp1(sorted_n, cdf_emp_sorted, cutoff, 'linear', 'extrap');
    
    current_share = interpolate - previous_cdf_value_emp;
    share_employment(i) = max(0, current_share);
    previous_cdf_value_emp = previous_cdf_value_emp + current_share;
end
share_employment(num_bins) = max(0, 1 - previous_cdf_value_emp);
share_employment = share_employment / sum(share_employment); % Re-normalize

% --- PLOT 4: FIRM SIZE PIE CHART ---
figure;
pie(share_firms, cellfun(@(x, y) sprintf('%s: %.1f%%', x, y*100), labels', num2cell(share_firms), 'UniformOutput', false));
title('Share of Firms by Number of Employees');

% --- PLOT 5: EMPLOYMENT SHARE PIE CHART ---
figure;
pie(share_employment, cellfun(@(x, y) sprintf('%s: %.1f%%', x, y*100), labels', num2cell(share_employment), 'UniformOutput', false));
title('Employment Share by Firm Size');