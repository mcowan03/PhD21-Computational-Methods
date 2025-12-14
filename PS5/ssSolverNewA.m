function [SS] = ssSolverNewA(A)

% Parameters
sigma = 1;
rho = 0.95;
sigma_epsilon = 0.009;
alpha = 0.33;
beta = 0.99;
delta = 0.025;
chi = 1;
mu = 0.60;
eta = 7.8827;

% At SS K'=K, a'=a, N=L, Psi=0

r_ss = 1/beta - 1;
K_N_ratio = ( (A * alpha) / (r_ss + delta) )^(1 / (1 - alpha));
w_ss = A * (1 - alpha) * K_N_ratio^alpha;
N_ss = sqrt(w_ss/(eta*(r_ss*K_N_ratio + w_ss)));
K_ss = K_N_ratio * N_ss; % Calculate steady state capital


Y_ss = A * K_ss^alpha * N_ss^(1-alpha);
I_ss = delta * K_ss;
c_ss = Y_ss - I_ss;

disp('--- Calibration Complete ---');
disp(['Steady State K: ', num2str(K_ss)]);
disp(['Steady State N: ', num2str(N_ss)]); 


% 5. Store in struct
    SS.K = K_ss; SS.N = N_ss; SS.Y = Y_ss;
    SS.C = c_ss; SS.I = I_ss; SS.w = w_ss; SS.r = r_ss;
end