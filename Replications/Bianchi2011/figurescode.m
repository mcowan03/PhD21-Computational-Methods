%% BIANCHI (2011) - CE vs SPP OVERLAY GRAPHICS
clear; clc; close all;

BURNIN = 500;
% 1. Load Data
ce  = load('prelimresults.mat');
spp = load('prelimresultsSPP.mat');

% Define the Crisis Shock Index (Lowest TFP/Endowment state)
crisis_idx = 1; 


[ce_b, ce_bp, ce_c, ce_css]   = extract_smooth_policy(ce, crisis_idx);
[spp_b, spp_bp, spp_c, spp_css] = extract_smooth_policy(spp, crisis_idx);

%% FIGURE 1: ERGODIC DISTRIBUTION OF DEBT
figure('Name', 'Ergodic Distribution of Debt', 'Units', 'inches', 'Position', [1, 1, 8, 5]);
hold on; grid on;

% Extract final valid paths
b_final_ce  = ce.vB(ce.BURNIN:end);
b_final_spp = spp.vB(spp.BURNIN:end);

% Kernel Density for smooth overlapping
[f_ce, xi_ce]   = ksdensity(b_final_ce);
[f_spp, xi_spp] = ksdensity(b_final_spp);

% Plot
fill(xi_ce, f_ce, [0.8 0.2 0.2], 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'LineWidth', 1.5);
fill(xi_spp, f_spp, [0.2 0.4 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'b', 'LineWidth', 1.5);

% Vertical lines for means
xline(mean(b_final_ce), 'r--', 'LineWidth', 2);
xline(mean(b_final_spp), 'b--', 'LineWidth', 2);

title('Ergodic Distribution of Debt Holdings', 'FontSize', 14);
xlabel('Debt Level (b)', 'FontSize', 12);
ylabel('Density', 'FontSize', 12);
legend('Competitive Equilibrium (CE)', 'Social Planner (SPP)', 'Mean CE', 'Mean SPP', 'Location', 'NorthWest');
hold off;

%% BACKING OUT POLICY FUNCTIONS
% The RTM solution is a "cloud" of points. We visualize it by slicing
% the cloud for specific shock realizations.

% 1. Setup Data
b_today = vB(BURNIN:end-1);
b_tomorrow = vB(BURNIN+1:end);
c_today = vC_T(BURNIN:end-1);
state_idx = vStatePath(BURNIN:end-1);

% 2. Select Representative Shocks to Plot (Dynamically adapts to Ny)
% Based on diagonal yT=yN
mid_val = ceil(Ny/2); % Finds the middle index dynamically

low_shock_idx  = 1;                                  % Lowest (yT=1, yN=1)
mid_shock_idx  = (mid_val-1)*Ny + mid_val;           % Median
high_shock_idx = Ny*Ny;                              % Highest (yT=Ny, yN=Ny)

figure('Name', 'Backing Out Policy Functions', 'Units', 'inches', 'Position', [1, 1, 12, 5]);

% --- PLOT BOND POLICY b'(b) ---
subplot(1,2,1); hold on; grid on;
title('Bond Policy Function b''(b, y)', 'FontSize', 12);
xlabel('Current Debt b_t', 'FontSize', 11); 
ylabel('Next Debt b_{t+1}', 'FontSize', 11);

% 45-degree line
plot([min(b_today), max(b_today)], [min(b_today), max(b_today)], 'k--', 'LineWidth', 1.5);

% Plot policies for specific shocks
plot_policy(b_today, b_tomorrow, state_idx, low_shock_idx, 'r', 'Low Shock');
plot_policy(b_today, b_tomorrow, state_idx, mid_shock_idx, 'g', 'Mid Shock');
plot_policy(b_today, b_tomorrow, state_idx, high_shock_idx, 'b', 'High Shock');
legend('45-degree', 'Low (Crisis)', 'Mid', 'High', 'Location', 'NorthWest');

% --- PLOT CONSUMPTION POLICY c(b) ---
subplot(1,2,2); hold on; grid on;
title('Consumption Policy Function c(b, y)', 'FontSize', 12);
xlabel('Current Debt b_t', 'FontSize', 11); 
ylabel('Consumption c^T_t', 'FontSize', 11);

plot_policy(b_today, c_today, state_idx, low_shock_idx, 'r', 'Low Shock');
plot_policy(b_today, c_today, state_idx, mid_shock_idx, 'g', 'Mid Shock');
plot_policy(b_today, c_today, state_idx, high_shock_idx, 'b', 'High Shock');
legend('Low (Crisis)', 'Mid', 'High', 'Location', 'NorthWest');




%% FIGURE 3: CONDITIONAL SADDLE PATHS (OVERLAY)
figure('Name', 'Conditional Saddle Paths', 'Units', 'inches', 'Position', [1, 1, 8, 6]);
hold on; grid on;

% CE Saddle Path
plot(ce_b, ce_c, 'r-', 'LineWidth', 2);
plot(ce_css.b, ce_css.c, 'pentagram', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');

% SPP Saddle Path
plot(spp_b, spp_c, 'b-', 'LineWidth', 2);
plot(spp_css.b, spp_css.c, 'pentagram', 'MarkerSize', 15, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');

% Formatting
title('Conditional Saddle Paths (Crisis State)', 'FontSize', 14);
xlabel('Current Debt (b_t)', 'FontSize', 12);
ylabel('Tradable Consumption (c^T)', 'FontSize', 12);
legend('CE Saddle Path', 'CE Steady State', 'SPP Saddle Path', 'SPP Steady State', 'Location', 'NorthWest');
hold off;

%% LOCAL HELPER FUNCTION
function [b_smooth, bp_smooth, c_smooth, css] = extract_smooth_policy(data, target_state)
    % Extracts and smooths the RTM scatter cloud for clear plotting
    valid = data.BURNIN:(data.pathlength-1);
    idx = valid(data.vStatePath(valid) == target_state);
    
    B_raw  = data.vB(idx);
    C_raw  = data.vC_T(idx);
    Bp_raw = data.vB(idx + 1);
    
    % Bin-average for smooth lines
    nBins = 100;
    B_lo = min(B_raw);  B_hi = max(B_raw);
    if B_hi - B_lo < 1e-6; B_lo = B_lo - 1e-4; B_hi = B_hi + 1e-4; end
    
    edges = linspace(B_lo, B_hi, nBins+1);
    centers = (edges(1:end-1) + edges(2:end)) / 2;
    C_bin = nan(nBins, 1); Bp_bin = nan(nBins, 1);
    
    for ib = 1:nBins
        in_bin = (B_raw >= edges(ib)) & (B_raw <= edges(ib+1));
        if sum(in_bin) >= 1
            C_bin(ib) = mean(C_raw(in_bin));
            Bp_bin(ib) = mean(Bp_raw(in_bin));
        end
    end
    
    keep = ~isnan(C_bin);
    b_smooth  = centers(keep)';
    c_smooth  = C_bin(keep);
    bp_smooth = Bp_bin(keep);
    
    % Calculate CSS
    B_curr = mean(B_raw);
    if length(b_smooth) > 1
        for ii = 1:5000
            B_curr = max(min(B_curr, B_hi), B_lo);
            B_next = interp1(b_smooth, bp_smooth, B_curr, 'linear', 'extrap');
            if abs(B_next - B_curr) < 1e-6; break; end
            B_curr = B_next;
        end
    end
    css.b = B_curr;
    css.c = interp1(b_smooth, c_smooth, B_curr, 'linear', 'extrap');
end

%% LOCAL HELPER FUNCTION
function plot_policy(x_data, y_data, states, target_state, color, name)
    % Filter the data for the specific shock
    idx = (states == target_state);
    x_target = x_data(idx);
    y_target = y_data(idx);
    
    % If there's no data for this state, exit gracefully
    if length(x_target) < 2
        return;
    end
    
    % Bin-average the scatter cloud to draw a smooth policy line
    nBins = 100;
    x_lo = min(x_target); 
    x_hi = max(x_target);
    
    % Safety catch for zero spread
    if x_hi - x_lo < 1e-6
        x_lo = x_lo - 1e-4; 
        x_hi = x_hi + 1e-4; 
    end
    
    edges = linspace(x_lo, x_hi, nBins+1);
    centers = (edges(1:end-1) + edges(2:end)) / 2;
    y_bin = nan(nBins, 1);
    
    for ib = 1:nBins
        in_bin = (x_target >= edges(ib)) & (x_target <= edges(ib+1));
        if sum(in_bin) >= 1
            y_bin(ib) = mean(y_target(in_bin));
        end
    end
    
    keep = ~isnan(y_bin);
    x_smooth = centers(keep)';
    y_smooth = y_bin(keep);
    
    % Plot the smoothed line
    plot(x_smooth, y_smooth, 'Color', color, 'LineWidth', 2.5, 'DisplayName', name);
end