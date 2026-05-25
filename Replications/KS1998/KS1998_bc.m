%% KS1998 Recursive Competitive Equilibrium (RCE) %%
close all; clc; clear variables;
tic;
%===================
% Load SS Solution %
%===================
load('ks1998_ss.mat');

ssK      = results.K;
ssV      = results.V;
ssPol    = results.pol_aprime;
ssDist   = results.dist;

%===================
% Parameters %
%===================
pAlpha = 0.36;
pBeta  = 0.99;
pDelta = 0.025;

%===================
% Individual Grids %
%===================
pnumgrida = 100;
vgridamin = 0;
vgridamax = 300;
x         = linspace(0, 0.5, pnumgrida);
y         = x.^5 / max(x.^5);
vgrida    = vgridamin + (vgridamax - vgridamin) * y;

pnumgridz = 2;
vgridz    = [0.25, 1.00];

%===================
% Transition Matrix — Full 4x4
%===================
% Rows/cols: (bad,unemp) (bad,emp) (good,unemp) (good,emp)
mtransz0 = [0.525,0.350,0.03125,0.09375;...
            0.035,0.84,0.0025,0.1225;...
            0.09375,0.03125,0.292,0.583;...
            0.0099,0.1151,0.0245,0.8505];

% Aggregate transition: mtransA(iA, iAp) = prob of going from A to A'
% Extracted from the 4x4 by summing blocks
mtransA_raw = [sum(sum(mtransz0(1:2,1:2))), sum(sum(mtransz0(1:2,3:4)));...
               sum(sum(mtransz0(3:4,1:2))), sum(sum(mtransz0(3:4,3:4)))];
mtransA(1,:) = mtransA_raw(1,:) / sum(mtransA_raw(1,:));
mtransA(2,:) = mtransA_raw(2,:) / sum(mtransA_raw(2,:));

% Idiosyncratic transition conditional on (iA -> iAp)
% Pi_zz{iA, iAp} is pnumgridz x pnumgridz
Pi_zz = cell(2,2);

% bad -> bad
tmp = mtransz0(1:2, 1:2);
Pi_zz{1,1}(1,:) = tmp(1,:) / sum(tmp(1,:));
Pi_zz{1,1}(2,:) = tmp(2,:) / sum(tmp(2,:));

% bad -> good
tmp = mtransz0(1:2, 3:4);
Pi_zz{1,2}(1,:) = tmp(1,:) / sum(tmp(1,:));
Pi_zz{1,2}(2,:) = tmp(2,:) / sum(tmp(2,:));

% good -> bad
tmp = mtransz0(3:4, 1:2);
Pi_zz{2,1}(1,:) = tmp(1,:) / sum(tmp(1,:));
Pi_zz{2,1}(2,:) = tmp(2,:) / sum(tmp(2,:));

% good -> good
tmp = mtransz0(3:4, 3:4);
Pi_zz{2,2}(1,:) = tmp(1,:) / sum(tmp(1,:));
Pi_zz{2,2}(2,:) = tmp(2,:) / sum(tmp(2,:));

%===================
% Aggregate Shock Grid %
%===================
pnumgridA = 2;
vgridA    = [0.99, 1.01];   % bad, good TFP

%===================
% Labour Supply by Aggregate State %
%===================
% Stationary distribution of z conditional on each aggregate state
% From eigenvectors of the conditional blocks
[vL_bad,~]  = eigs(Pi_zz{1,1}', 1);
vL_bad      = real(vL_bad) / sum(real(vL_bad));
[vL_good,~] = eigs(Pi_zz{2,2}', 1);
vL_good     = real(vL_good) / sum(real(vL_good));

supplyL     = zeros(1, pnumgridA);
supplyL(1)  = vgridz * vL_bad;    % labour in bad state
supplyL(2)  = vgridz * vL_good;   % labour in good state

fprintf('supplyL: bad=%.4f good=%.4f\n', supplyL(1), supplyL(2));

%===================
% Aggregate Capital Grid %
%===================
pnumgridK = 8;
K_dev     = 0.25;
vgridK    = linspace(ssK*(1-K_dev), ssK*(1+K_dev), pnumgridK);

%===================
% Simulation Settings %
%===================
T_sim      = 3000;
T_burn     = 500;
rng(42);

% Draw aggregate shock sequence
A_idx    = zeros(T_sim, 1);
A_idx(1) = 1;   % start in bad state
for t = 2:T_sim
    if rand < mtransA(A_idx(t-1), 1)
        A_idx(t) = 1;
    else
        A_idx(t) = 2;
    end
end

%===================
% LOM Initial Guess %
%===================
lom_a = [0; 0];   % [a_bad; a_good]
lom_b = [1; 1];   % [b_bad; b_good]

%===================
% VFI Initialisation — Warm Start from SS %
%===================
V          = zeros(pnumgrida, pnumgridz, pnumgridK, pnumgridA);
pol_aprime = zeros(pnumgrida, pnumgridz, pnumgridK, pnumgridA);
for iA = 1:pnumgridA
    for iK = 1:pnumgridK
        V(:,:,iK,iA)          = ssV;
        pol_aprime(:,:,iK,iA) = ssPol;
    end
end

%===================
% Numerical Settings %
%===================
tol_lom      = 1e-4;
tol_vfi      = 1e-6;
weight_lom   = 0.5;
max_lom_iter = 100;

%===================
% OUTER LOOP: LOM %
%===================
error_lom = 10;
iter_lom  = 0;

fprintf('Starting KS1998 dynamic solution...\n\n');

while error_lom > tol_lom && iter_lom < max_lom_iter
    iter_lom  = iter_lom + 1;
    lom_a_old = lom_a;
    lom_b_old = lom_b;
    fprintf('=== LOM iteration %d ===\n', iter_lom);

    %===================
    % STEP 1: VFI
    %===================
    errorvfi = 10;
    iter_vfi = 0;

    while errorvfi > tol_vfi && iter_vfi < 2000
        iter_vfi = iter_vfi + 1;
        V_prev   = V;

        for iA = 1:pnumgridA
            A = vgridA(iA);
            L = supplyL(iA);

            for iK = 1:pnumgridK
                Kbar = vgridK(iK);

                % Prices
                r = pAlpha * A * (Kbar/L)^(pAlpha-1) - pDelta;
                w = (1-pAlpha) * A * (Kbar/L)^pAlpha;

                % Next period Kbar from LOM
                logKbar_next = lom_a(iA) + lom_b(iA) * log(Kbar);
                Kbar_next    = exp(logKbar_next);
                Kbar_next    = min(max(Kbar_next, vgridK(1)), vgridK(end));

                % Expected continuation value EV(ia, iz)
                EV = zeros(pnumgrida, pnumgridz);
                for iAp = 1:pnumgridA
                    for iz = 1:pnumgridz
                        for izp = 1:pnumgridz

                            % Interpolate V(a, izp, Kbar_next, iAp) over K grid
                            V_at_Kp = interp1(vgridK, ...
                                squeeze(V(:, izp, :, iAp))', ...
                                Kbar_next, 'linear', 'extrap')';
                            % pnumgrida x 1

                            % Weight by idiosyncratic and aggregate transitions
                            EV(:, iz) = EV(:, iz) + ...
                                mtransA(iA, iAp) * Pi_zz{iA,iAp}(iz, izp) * V_at_Kp;
                        end
                    end
                end

                % Discounted continuation value
                CV = pBeta * EV;   % pnumgrida x pnumgridz

                % Maximise for each (ia, iz)
                for iz = 1:pnumgridz
                    CV_col = CV(:, iz);

                    for ia = 1:pnumgrida
                        coh = (1+r)*vgrida(ia) + w*vgridz(iz);

                        % Grid search
                        V_best  = -1e10;
                        ap_best = 1;
                        for iap = 1:pnumgrida
                            c = coh - vgrida(iap);
                            if c <= 0, continue; end
                            vj = log(c) + CV_col(iap);
                            if vj > V_best
                                V_best  = vj;
                                ap_best = iap;
                            end
                        end

                        % fminbnd refinement
                        ilo  = max(ap_best-1, 1);
                        ihi  = min(ap_best+1, pnumgrida);
                        a_lo = vgrida(ilo);
                        a_hi = min(vgrida(ihi), coh - 1e-8);

                        if a_hi > a_lo + 1e-10
                            obj = @(ap) -(log(coh-ap) + interp1(vgrida, CV_col, ap, 'linear'));
                            [apopt, negval]         = fminbnd(obj, a_lo, a_hi);
                            pol_aprime(ia,iz,iK,iA) = apopt;
                            V(ia,iz,iK,iA)          = -negval;
                        else
                            pol_aprime(ia,iz,iK,iA) = vgrida(ap_best);
                            V(ia,iz,iK,iA)          = V_best;
                        end
                    end
                end
            end
        end

        errorvfi = max(abs(V(:) - V_prev(:)));
        if mod(iter_vfi, 50) == 0
            fprintf('  VFI iter %d | error = %.2e\n', iter_vfi, errorvfi);
        end
    end
    fprintf('  VFI converged in %d iters\n', iter_vfi);

    %===================
    % STEP 2: Simulate
    %===================
    fprintf('  Simulating %d periods...\n', T_sim);

    Kbar_sim    = zeros(T_sim, 1);
    Kbar_sim(1) = ssK;
    currentdist = ssDist;

    for t = 1:T_sim-1
        iA   = A_idx(t);
        iAp  = A_idx(t+1);
        Kbar = Kbar_sim(t);

        nextdist = zeros(pnumgrida, pnumgridz);

        for iz = 1:pnumgridz
            for ia = 1:pnumgrida
                mass = currentdist(ia, iz);
                if mass < 1e-12, continue; end

                % Policy at this Kbar (interpolate over K grid)
                nexta = interp1(vgridK, squeeze(pol_aprime(ia,iz,:,iA)), ...
                    Kbar, 'linear', 'extrap');
                nexta = min(max(nexta, vgrida(1)), vgrida(end));

                % Bracket for histogram
                lb = sum(vgrida < nexta);
                lb = max(lb, 1);
                lb = min(lb, pnumgrida-1);
                ub = lb + 1;

                wlb = (vgrida(ub) - nexta) / (vgrida(ub) - vgrida(lb));
                wlb = min(max(wlb, 0), 1);
                wub = 1 - wlb;

                % Use conditional idiosyncratic transition Pi_zz{iA, iAp}
                for izp = 1:pnumgridz
                    prob = Pi_zz{iA, iAp}(iz, izp);
                    nextdist(lb, izp) = nextdist(lb, izp) + mass * prob * wlb;
                    nextdist(ub, izp) = nextdist(ub, izp) + mass * prob * wub;
                end
            end
        end

        currentdist  = nextdist;
        Kbar_sim(t+1) = vgrida(:)' * sum(currentdist, 2);
    end

    %===================
    % STEP 3: OLS Regression to update LOM
    %===================
    Kbar_t  = Kbar_sim(T_burn:T_sim-1);
    Kbar_tp = Kbar_sim(T_burn+1:T_sim);
    A_t     = A_idx(T_burn:T_sim-1);

    lom_a_new = zeros(pnumgridA, 1);
    lom_b_new = zeros(pnumgridA, 1);
    R2        = zeros(pnumgridA, 1);

    for iA = 1:pnumgridA
        idx  = (A_t == iA);
        logK  = log(Kbar_t(idx));
        logKp = log(Kbar_tp(idx));

        X        = [ones(sum(idx),1), logK];
        beta_ols = X \ logKp;

        lom_a_new(iA) = beta_ols(1);
        lom_b_new(iA) = beta_ols(2);

        logKp_hat = X * beta_ols;
        SS_res    = sum((logKp - logKp_hat).^2);
        SS_tot    = sum((logKp - mean(logKp)).^2);
        R2(iA)    = 1 - SS_res / SS_tot;
    end

    fprintf('  Bad:  a=%.6f  b=%.6f  R2=%.8f\n', lom_a_new(1), lom_b_new(1), R2(1));
    fprintf('  Good: a=%.6f  b=%.6f  R2=%.8f\n', lom_a_new(2), lom_b_new(2), R2(2));

    %===================
    % STEP 4: Update LOM
    %===================
    error_lom = max(abs([lom_a_new - lom_a_old; lom_b_new - lom_b_old]));
    fprintf('  LOM error = %.2e\n\n', error_lom);

    lom_a = weight_lom * lom_a_old + (1-weight_lom) * lom_a_new;
    lom_b = weight_lom * lom_b_old + (1-weight_lom) * lom_b_new;
end

%===================
% Report %
%===================
if error_lom < tol_lom
    fprintf('\nKS1998 Dynamic Solved!\n');
else
    fprintf('\nDid not converge after %d iterations\n', max_lom_iter);
end
toc;
fprintf('LOM bad:  log(K'') = %.6f + %.6f*log(K)\n', lom_a(1), lom_b(1));
fprintf('LOM good: log(K'') = %.6f + %.6f*log(K)\n', lom_a(2), lom_b(2));

%===================
% Store Results %
%===================
dynresults.lom_a      = lom_a;
dynresults.lom_b      = lom_b;
dynresults.R2         = R2;
dynresults.V          = V;
dynresults.pol_aprime = pol_aprime;
dynresults.Kbar_sim   = Kbar_sim;
dynresults.A_idx      = A_idx;
dynresults.vgridK     = vgridK;
dynresults.vgridA     = vgridA;
dynresults.vgrida     = vgrida;
dynresults.vgridz     = vgridz;
dynresults.Pi_zz      = Pi_zz;
dynresults.mtransA    = mtransA;
dynresults.supplyL    = supplyL;

save('ks1998_dynamic.mat', 'dynresults');
fprintf('Saved to ks1998_dynamic.mat\n');

%===================
% Plots %
%===================
figure;
subplot(2,1,1);
plot(Kbar_sim);
xlabel('Period'); ylabel('K');
title('Simulated Aggregate Capital');
grid on;

subplot(2,1,2);
plot(A_idx);
xlabel('Period'); ylabel('A (1=bad, 2=good)');
title('Aggregate Shock Sequence');
grid on;