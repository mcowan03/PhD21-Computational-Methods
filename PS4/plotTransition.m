% --- TCE PLOTS ---

% Note: I run all TCE codes and store K_t paths with different names.
figure;
plot(1:T, K_path_decay, 'b-', 'LineWidth', 2);
hold on; % Call this once after the first plot
plot(1:T, K_path_transitive, 'g-', 'LineWidth', 2);
plot(1:T, K_path_newSS, 'r-', 'LineWidth', 2);
plot([1, T], [K_old, K_old], 'k--', 'LineWidth', 1);
plot([1, T], [K_new, K_new], 'm--', 'LineWidth', 1);
title('Aiyagari Perfect Foresight Transition');
xlabel('Time Period');
ylabel('Aggregate Capital (K_t)');


% Provide 5 labels, one for each of the 5 plot commands above
legend('Decaying Shock Path', ...
       'Transitory Shock Path', ... 
       'Permanent Shock Path', ... 
       sprintf('Old SS K*=%.4f', K_old), ...
       sprintf('New SS K*=%.4f', K_new), ...
       'Location', 'best'); 

grid on;
box on;
set(gca, 'FontSize', 12);
fprintf('\nTransition solved! Final error: %.4e\n', error_path);