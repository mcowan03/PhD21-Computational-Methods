%% Finding eta & SS values%%

function [SS, eta] = ssSolver(A)

% Parameters
sigma = 1;
rho = 0.95;
sigma_epsilon = 0.009;
alpha = 0.33;
beta = 0.99;
delta = 0.025;
chi = 1;
mu = 0.60;



% N target
N_target = 1/3;

% At SS K'=K, a'=a, N=L, Psi=0

r_ss = 1/beta - 1;
K_N_ratio = ( (A * alpha) / (r_ss + delta) )^(1 / (1 - alpha));

N_ss = N_target;
K_ss = K_N_ratio * N_ss;
w_ss = A * (1 - alpha) * K_N_ratio^alpha;
Y_ss = A * K_ss^alpha * N_ss^(1-alpha);
I_ss = delta * K_ss;
c_ss = Y_ss - I_ss;

eta_calibrated = w_ss/(N_ss*c_ss);
eta = eta_calibrated;

disp('--- Calibration Complete ---');
disp(['Calibrated eta: ', num2str(eta_calibrated)]);
disp(['Steady State K: ', num2str(K_ss)]);
disp(['Steady State N: ', num2str(N_ss)]); 


% 5. Store in struct
    SS.K = K_ss; SS.N = N_ss; SS.Y = Y_ss;
    SS.C = c_ss; SS.I = I_ss; SS.w = w_ss; SS.r = r_ss;
end