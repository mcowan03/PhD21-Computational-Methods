%% 7. Data Visualization: Life-Cycle Plots

% --- 1. Compute Averages and Standard Deviations ---

% The Sim_ arrays have dimensions (Time, Household)
% The mean across the household dimension (dim=2) gives the average life-cycle path.
Mean_A = mean(Sim_A, 2, 'omitnan');
Mean_H = mean(Sim_H, 2, 'omitnan');
Mean_C = mean(Sim_C, 2, 'omitnan');
Mean_L = mean(Sim_L, 2, 'omitnan'); % This is the Labor Force Participation Rate

% Calculate standard deviation of Assets to show wealth inequality
Std_A = std(Sim_A, 0, 2, 'omitnan');

% Time vector for the x-axis (t=0 to t=T)
Time_Vec = 0:T;

% --- 2. Create the Figure and Subplots ---

figure('Name', 'Life-Cycle Simulation Results', 'Position', [100, 100, 1000, 800]);

% --- Subplot 1: Assets (A) and Wealth Inequality ---
subplot(2, 2, 1);
plot(Time_Vec, Mean_A, 'LineWidth', 2, 'DisplayName', 'Average Assets');
hold on;
plot(Time_Vec, Mean_A + 2*Std_A, '--', 'Color', [0.5 0.5 0.5], 'DisplayName', 'Mean + 2 Std. Dev.');
plot(Time_Vec, Mean_A - 2*Std_A, '--', 'Color', [0.5 0.5 0.5], 'DisplayName', 'Mean - 2 Std. Dev.');
hold off;
title(' Assets and Wealth Distribution');
xlabel('Age (t)');
ylabel('Asset Level');
grid on;
legend('show', 'Location', 'northwest');


% --- Subplot 2: Human Capital (H) ---
subplot(2, 2, 2);
plot(Time_Vec, Mean_H, 'LineWidth', 2, 'Color', 'g');
title(' Human Capital Accumulation');
xlabel('Age (t)');
ylabel('Human Capital Level (h)');
grid on;


% --- Subplot 3: Consumption (C) ---
subplot(2, 2, 3);
plot(Time_Vec, Mean_C, 'LineWidth', 2, 'Color', 'r');
title(' Consumption Path');
xlabel('Age (t)');
ylabel('Consumption (c)');
grid on;


% --- Subplot 4: Labor Force Participation (L) ---
subplot(2, 2, 4);
plot(Time_Vec, Mean_L, 'LineWidth', 2, 'Color', 'b');
title(' Labor Force Participation Rate');
xlabel('Age (t)');
ylabel('Probability of Working');
ylim([0, 1]); % Probability must be between 0 and 1
grid on;

% --- Final Touches ---
sgtitle(['Simulated Life-Cycle Paths (N=', num2str(N_households), ')']);

