%% Bianchi (2011) - Social Planner Problem (SPP) %%
clear; clc; close all;

% Parameters
r = 0.04;
sigma = 2;
eta = 1/0.83 - 1;
kappaN = 0.32;
kappaT = 0.32;
omega = 0.31;
beta = 0.91;

%% EXOGENOUS SHOCKS %%
Ny = 4;
rho_y = 0.8;
sigma_y = 0.058;
y_bar = 0;
nStd = 3;

[logygrid, P] = tauchen(Ny, rho_y, sigma_y, y_bar, nStd);
yTGrid = exp(logygrid);
yNGrid = exp(logygrid);

% Simulation
T = 5000;
BURNIN = 0.1*T;
pathlength = T + BURNIN;
maxiter = 15000;
tol = 1e-5;
damp = 0.99;

%% SIMULATING A PATH FOR EXOGENOUS SHOCKS
iniPoint = 1;
yTSimPath = fnSimulator(iniPoint,P,pathlength);
yNSimPath = fnSimulator(iniPoint,P,pathlength);

%% Initial guess for allocation path
b0 = -0.2;
cTss = 1 + b0*r;

vB = b0*ones(pathlength, 1) + normrnd(0,0.0000001,pathlength,1); 
vC_T = cTss*ones(pathlength, 1);    
vMu = zeros(pathlength, 1); 

% SPP SPECIFIC: Initialize Planner's Lambda
vYN = yNGrid(yNSimPath)';
vYT = yTGrid(yTSimPath)';
vLambda_sp = get_uT_deriv_val(vC_T, vYN, omega, eta, sigma); 

% Separate paths to be updated iteratively
vBnew = zeros(pathlength, 1);
vC_Tnew = zeros(pathlength, 1);
vMunew = zeros(pathlength, 1);
vLambda_sp_new = zeros(pathlength, 1);

%% RTM LOOP %%
tic;
pNumIter = 0;
error = 2;

% Pre calculate states
N_states = Ny*Ny;
vStatePath = (yTSimPath - 1) * Ny + yNSimPath; 
P_State = kron(P, P);
iS = vStatePath;
iFuture = [(2:pathlength)'; pathlength];
futureShock = vStatePath(iFuture);

mTransRealized = zeros(pathlength, 1);
for t = 1:pathlength
    mTransRealized(t) = P_State(iS(t), futureShock(t));
end

while error > 1e-4
    % Backwards solution
    vBprime = [vB(2:end); vB(1)];
    tempExpMU = zeros(pathlength, 1);
    
    for iStatePrime = 1:N_states
        candidateLoc = find(vStatePath == iStatePrime);
        candidateLoc(candidateLoc > pathlength - BURNIN) = [];
        candidateLoc(candidateLoc < BURNIN) = [];
        if length(candidateLoc) < 2; continue; end 
        
        candidate = vB(candidateLoc);
        [candidate, index] = sort(candidate);
        candidateLoc = candidateLoc(index);
        BLow = sum(vBprime > candidate', 2); 
        
        BLow(BLow <= 1) = 1;
        BLow(BLow >= length(index)) = length(index) - 1;
        BHigh = BLow + 1;
        
        candHigh = candidate(BHigh);
        candLow  = candidate(BLow);
        denom = candHigh - candLow;
        denom(denom == 0) = 1e-10; 
        weightLow = (candHigh - vBprime) ./ denom;
        weightLow(weightLow < 0) = 0;
        weightLow(weightLow > 1) = 1;
        
        % SPP SPECIFIC: Interpolate Planner's Lambda, not just u'(c)
        MU_Low  = vLambda_sp(candidateLoc(BLow));
        MU_High = vLambda_sp(candidateLoc(BHigh));
        MU_Interp = weightLow .* MU_Low + (1 - weightLow) .* MU_High;
        
        prob_trans = P_State(iS, iStatePrime);
        mask = (futureShock ~= iStatePrime); 
        tempExpMU = tempExpMU + mask .* prob_trans .* MU_Interp;
    end
    
    MU_Realized = vLambda_sp(iFuture); 
    tempExpMU = tempExpMU + mTransRealized .* MU_Realized;
    
    % SPP EULER RHS: beta * (1+r) * E_t[\lambda^{sp}_{t+1}]
    rhs_euler = beta * (1+r) * tempExpMU;
    
    % --- STEP 2: FORWARD SIMULATION ---
    for t = 1:pathlength
        yt = vYT(t);
        yn = vYN(t);
        bt = vB(t);
        rhs = rhs_euler(t);
        
        c_min = 1e-5;
        c_max = yt + (1+r)*max(abs(vB)) + 10; 
        
        % A. Unconstrained Solution (mu = 0, so lambda_sp = u'(c))
        resid_unc = @(c) get_uT_deriv_val(c, yn, omega, eta, sigma) - rhs;
        try
            if resid_unc(c_min) * resid_unc(c_max) > 0
                c_unc = yt; 
            else
                c_unc = fzero(resid_unc, [c_min, c_max]);
            end
        catch
            c_unc = yt; 
        end
        
        b_prime_unc = (1+r)*bt + yt - c_unc;
        
        % B. Check Constraint
        pN_unc = ((1-omega)/omega) * (c_unc/yn)^(eta+1);
        limit_val = -(kappaN * pN_unc * yn + kappaT * yt);
        
        if b_prime_unc >= limit_val
            % Not Binding
            vC_Tnew(t) = c_unc;
            vBnew(t+1) = b_prime_unc;
            vMunew(t)  = 0;
            vLambda_sp_new(t) = get_uT_deriv_val(c_unc, yn, omega, eta, sigma); % mu=0
        else
            % C. Binding Constraint
            resid_bind = @(c) ((1+r)*bt + yt - c) + ...
                              (kappaN * ((1-omega)/omega)*(c/yn)^(eta+1)*yn + kappaT * yt);
            try
                if resid_bind(1e-6) * resid_bind(c_unc) > 0
                     c_bind = fminbnd(@(c) abs(resid_bind(c)), 1e-6, c_unc);
                else
                    c_bind = fzero(resid_bind, [1e-6, c_unc]);
                end
            catch
                c_bind = c_unc; 
            end
            
            vC_Tnew(t) = c_bind;
            vBnew(t+1) = (1+r)*bt + yt - c_bind;
            
            % SPP SPECIFIC: Calculate mu based on Planner's Euler
            Psi_val = kappaN * ((1-omega)/omega) * (1+eta) * (c_bind/yn)^eta;
            u_prime_c = get_uT_deriv_val(c_bind, yn, omega, eta, sigma);
            
            % Derived from: u'(c) + mu*Psi = rhs + mu
            mu_val = (u_prime_c - rhs) / (1 - Psi_val);
            vMunew(t) = max(0, mu_val);
            
            % Update Planner's Lambda
            vLambda_sp_new(t) = u_prime_c + vMunew(t) * Psi_val;
        end
    end
    
    vBnew(1) = vB(1); 
    vBnew = vBnew(1:pathlength);
    
    % --- STEP 3: UPDATE ---
    diff = max(abs(vB - vBnew));
    error = diff;
    
    vB         = damp*vB   + (1-damp)*vBnew;
    vC_T       = damp*vC_T + (1-damp)*vC_Tnew;
    vMu        = damp*vMu  + (1-damp)*vMunew;
    vLambda_sp = damp*vLambda_sp + (1-damp)*vLambda_sp_new; % Update SPP Lambda
    % After updating vC_T and vMu, recompute lambda analytically
    % Safeguard before recomputing lambda
%vC_T = max(vC_T, 1e-8); %NEW
% vYN  = max(vYN,  1e-8); % NEW
   % Psi_path = kappaN * ((1-omega)/omega) * (1+eta) * (vC_T ./ vYN).^eta; % NEW
   % vLambda_sp = get_uT_deriv_val(vC_T, vYN, omega, eta, sigma) + vMu .* Psi_path; % NEW

    if (floor((pNumIter-1)/50) == (pNumIter-1)/50)
        fprintf('Iter: %d | Error: %.6f\n', pNumIter, error);
    end
    
    pNumIter = pNumIter + 1;
    if pNumIter > maxiter; break; end
end
toc;

%% STATISTICS %%

%% CALCULATE CRISIS FREQUENCY
% Remove the burn-in period to calculate the true ergodic probability
mu_final = vMu(BURNIN:end);

% Define a small threshold to filter out numerical solver noise 
% (e.g., a solver might return mu = 1e-12 instead of exactly 0)
crisis_threshold = 1e-5; 

% Identify periods where the constraint binds
crisis_indicator = (mu_final > crisis_threshold);

% Calculate the frequency
crisis_frequency = sum(crisis_indicator) / length(mu_final);

fprintf('\n==================================================\n');
fprintf('Model Diagnostics:\n');
fprintf('Crisis Frequency: %.2f%%\n', crisis_frequency * 100);
fprintf('==================================================\n');


%% Consumption %%
%% CALCULATE AVERAGE CONSUMPTION
% Remove the burn-in period to calculate true ergodic moments
cT_final = vC_T(BURNIN:end);
cN_final = vYN(BURNIN:end)'; % Market clearing: cN = yN

% 1. Average Tradable Consumption
mean_cT = mean(cT_final);

% 2. Average Non-Tradable Consumption 
mean_cN = mean(cN_final);

% 3. Average Aggregate Consumption (CES Basket)
% Formula: c = [omega * (cT)^(-eta) + (1-omega) * (cN)^(-eta)]^(-1/eta)
c_basket_final = (omega * cT_final.^(-eta) + (1-omega) * cN_final.^(-eta)).^(-1/eta);
mean_c_basket = mean(c_basket_final);

%% CALCULATE CONSUMPTION VOLATILITY
% Calculate Standard Deviation (Volatility)
std_cT = std(cT_final);
std_cN = std(cN_final);
std_c_basket = std(c_basket_final);

%% CALCULATE AVERAGE DEBT
% Remove the burn-in period to calculate the true ergodic mean
b_final = vB(BURNIN:end);

% Calculate Average Debt
mean_b = mean(b_final);

% Print the result
fprintf('\n==================================================\n');
fprintf('Debt Diagnostics:\n');
fprintf('Average Debt (b): %.4f\n', mean_b);
fprintf('==================================================\n');

% Print the results
fprintf('\n==================================================\n');
fprintf('SOCIAL PLANNER ! Consumption Volatility (Standard Deviation):\n');
fprintf('Volatility of Tradable Consumption (cT): %.4f\n', std_cT);
fprintf('Volatility of Non-Tradable Consumption (cN): %.4f\n', std_cN);
fprintf('Volatility of Aggregate Consumption (c): %.4f\n', std_c_basket);
fprintf('==================================================\n');
fprintf('\n==================================================\n');
fprintf('Consumption Diagnostics:\n');
fprintf('Average Tradable Consumption (cT): %.4f\n', mean_cT);
fprintf('Average Non-Tradable Consumption (cN): %.4f\n', mean_cN);
fprintf('Average Aggregate Consumption (c): %.4f\n', mean_c_basket);
fprintf('==================================================\n');

save('..\Bianchi2011\prelimresultsSPP.mat')

%% LOCAL HELPER FUNCTIONS %%
function [Z,Zprob] = tauchen(N,rho,sigma,mu,n_std)
    Z     = zeros(N,1);
    Zprob = zeros(N,N);
    a     = (1-rho)*mu;
    Z(N)  = n_std * sqrt(sigma^2 / (1 - rho^2));
    Z(1)  = -Z(N);
    step  = (Z(N) - Z(1)) / (N - 1);
    for i=2:(N-1)
        Z(i) = Z(1) + step * (i - 1);
    end 
    Z = Z + a / (1-rho);
    for j = 1:N
        for k = 1:N
            if k == 1
                Zprob(j,k) = normcdf((Z(1) - a - rho * Z(j) + step / 2) / sigma);
            elseif k == N
                Zprob(j,k) = 1 - normcdf((Z(N) - a - rho * Z(j) - step / 2) / sigma);
            else
                Zprob(j,k) = normcdf((Z(k) - a - rho * Z(j) + step / 2) / sigma) - ...
                             normcdf((Z(k) - a - rho * Z(j) - step / 2) / sigma);
            end
        end
    end
end

function simPath = fnSimulator(ini, P, T)
    simPath = zeros(T,1);
    simPath(1) = ini;
    cdf = cumsum(P,2);
    for t=2:T
        r = rand;
        simPath(t) = find(cdf(simPath(t-1),:) >= r, 1, 'first');
    end
end

function u_prime = get_uT_deriv_val(cT, cN, omega, eta, sigma)
    % Marginal utility of tradables using CES aggregator
    c = (omega*cT.^(-eta) + (1-omega)*cN.^(-eta)).^(-1/eta);
    dc_dcT = omega * (c ./ cT).^(eta+1);
    u_prime = (c.^(-sigma)) .* dc_dcT;
end