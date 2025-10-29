% Set the plotting range and resolution
w_min = 0.2;
w_max = 2;
num_points = 20;
w_grid = linspace(w_min, w_max, num_points);

% Set a fixed transfer T (e.g., the equilibrium value you found, or T=0)
T_fixed = 0.3757; % Example: Using your initial guess or final GE value

% Initialize vector to store results
L_agg_vector = zeros(size(w_grid));

% Run the simulation for each wage
for j = 1:num_points
    w_j = w_grid(j);
    L_agg_vector(j) = Calculate_Lagg(w_j, T_fixed); 
end

% Plot the aggregate labor supply curve
figure;
plot(L_agg_vector, w_grid, 'b-', 'LineWidth', 2);
xlabel('Aggregate Labor Supply (L_{agg})');
ylabel('Wage (w)');
title(['Aggregate Labor Supply Curve (Fixed T = ', num2str(T_fixed), ')']);
grid on;

hold on;
% Overlay labour demand...

alpha = 0.3;
A = 1;

labDemand = (((1-alpha)*A)./w_grid).^(1/alpha)

% Overlay labor demand curve
plot(labDemand, w_grid, 'r--', 'LineWidth', 2, 'DisplayName', 'Labor Demand');

% Add legend
legend('Labor Supply', 'Labor Demand');

% Release the hold
hold off;