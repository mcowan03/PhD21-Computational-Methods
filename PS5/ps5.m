%% ASSIGNMENT 5 %%

% Part 1:
% Turn off aggregate uncertainty, A=1;
% Solve for steady state and calibrate eta to match equilibrium labor hour
% = 1/3

% Parameters
sigma = 1;
rho = 0.95;
sigma_epsilon = 0.009;
alpha = 0.33;
beta = 0.99;
delta = 0.025;
chi = 1;
mu = 0.60;

%% Setup Grids %%
% Wealth grid (little k, or 'a')
a_lower = 0.0;
a_upper = 150.0;
Na = 150;

% Maliar, Maliar and Valli, 2010 finer grids at low a:
x = linspace(0,0.5,Na);
y = x.^5/max(x.^5);
a_grid = a_lower+(a_upper-a_lower)*y;
K_grid = a_grid;

%% Productivity %%
% No idiosyncratic productivity !

% Aggregate TFP
A = 1; % = 1 for now... We may be able to solve for steady state analytically

%% Solving for equilibrium %%
% Note that when A=1 the system has no uncertainty.
% Also note that the representative household does not internalise that
% a' = K'

% Pre-allocation
valFunc = zeros(Na, Na); % V(a; K')
aprimePol = zeros(Na, Na);
nPol = zeros(Na, Na);



adjustmentCost = (mu/2)*()^2*a

