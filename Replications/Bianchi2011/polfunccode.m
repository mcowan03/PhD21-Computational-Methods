%% BACKING OUT POLICY FUNCTIONS: CE vs SPP (MID SHOCK ONLY)
clear; clc; close all;

% 1. Load Data from both simulations
ce  = load('prelimresults.mat');
spp = load('prelimresultsSPP.mat');

% Extract shared parameters
BURNIN = ce.BURNIN;
Ny = ce.Ny;

% Setup CE Data
ce_b_today = ce.vB(BURNIN:end-1);
ce_b_tomorrow = ce.vB(BURNIN+1:end);
ce_c_today = ce.vC_T(BURNIN:end-1);
ce_state = ce.vStatePath(BURNIN:end-1);

% Setup SPP Data
spp_b_today = spp.vB(BURNIN:end-1);
spp_b_tomorrow = spp.vB(BURNIN+1:end);
spp_c_today = spp.vC_T(BURNIN:end-1);
spp_state = spp.vStatePath(BURNIN:end-1);

% 2. Select Representative Shock to Plot
% Mid Shock: Dynamically finding the median state
mid_val = ceil(Ny/2); 
mid_shock_idx = (mid_val-1)*Ny + mid_val;                              
mid_shock_idx = 6
figure('Name', 'Policy Functions: CE vs SPP (Mid Shock)', 'Units', 'inches', 'Position', [1, 1, 14, 6]);

% --- PLOT BOND POLICY b'(b) ---
subplot(1,2,1); hold on; grid on;
title('Bond Policy Function b''(b, y) - Mid Shock', 'FontSize', 12);
xlabel('Current Debt b_t', 'FontSize', 11); 
ylabel('Next Debt b_{t+1}', 'FontSize', 11);

% 45-degree line
plot([min(ce_b_today), max(ce_b_today)], [min(ce_b_today), max(ce_b_today)], 'k:', 'LineWidth', 1.5, 'DisplayName', '45-degree');

% CE Policy (Solid Line)
plot_policy(ce_b_today, ce_b_tomorrow, ce_state, mid_shock_idx, 'g', '-', 'CE - Mid Shock');

% SPP Policy (Dashed Line)
plot_policy(spp_b_today, spp_b_tomorrow, spp_state, mid_shock_idx, 'r', '--', 'SPP - Mid Shock');

legend('Location', 'NorthWest');

% --- PLOT CONSUMPTION POLICY c(b) ---
subplot(1,2,2); hold on; grid on;
title('Consumption Policy Function c(b, y) - Mid Shock', 'FontSize', 12);
xlabel('Current Debt b_t', 'FontSize', 11); 
ylabel('Consumption c^T_t', 'FontSize', 11);

% CE Policy (Solid Line)
plot_policy(ce_b_today, ce_c_today, ce_state, mid_shock_idx, 'g', '-', 'CE - Mid Shock');

% SPP Policy (Dashed Line)
plot_policy(spp_b_today, spp_c_today, spp_state, mid_shock_idx, 'r', '--', 'SPP - Mid Shock');

legend('Location', 'NorthWest');


%% LOCAL HELPER FUNCTION
function plot_policy(x_data, y_data, states, target_state, color, linestyle, name)
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
    
    % Plot the smoothed line with Dynamic LineStyle and Name mapping
    plot(x_smooth, y_smooth, 'Color', color, 'LineStyle', linestyle, 'LineWidth', 2.5, 'DisplayName', name);
end