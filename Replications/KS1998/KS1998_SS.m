%% KS1998 Stationary Recursive Competitive Equilibrium (SRCE) %%
close all;
clc;
clear variables;

%============
% Algorithm %
%============

%=============================================
    % 1. Compute stationary z distribution
    % 2. Compute agg. L
    % 3. Guess agg. K -> {r,w}
    % 4. Solve individual problem by VFI
    % 5. Compute implied K:
    %   5.1. Solve for invariant distribution  
    %       5.2. impK = int{a'(a,z,S)}dPhi
    % 6. Check error: abs(K - impK)
    % 7. Update and back to 3.
%=============================================



%===================
% Setup Parameters %
%===================

% Model Parameters 
pAlpha = 0.36;
pBeta = 0.99;
pDelta = 0.025;


% Numerical Parameters
weight = 0.8;


%idiosyncratic income shock
pnumgridz   = 2;
vgridz      = [0.25,1.00];
pnumgridA   = 2;
mtransz     = zeros(pnumgridz,pnumgridz);
mtransz0    = [0.525,0.350,0.03125,0.09375;...
               0.035,0.84,0.0025,0.1225;...
               0.09375,0.03125,0.292,0.583;...
               0.0099,0.1151,0.0245,0.8505];    
mtransA = [sum(sum(mtransz0(1:2,1:2))),sum(sum(mtransz0(1:2,3:4))); ...
           sum(sum(mtransz0(3:4,1:2))),sum(sum(mtransz0(3:4,3:4)))];
mtransA = mtransA/(2*pnumgridA);
mtransz0(1,1:2) = mtransz0(1,1:2)/sum(mtransz0(1,1:2),'all');
mtransz0(2,1:2) = mtransz0(2,1:2)/sum(mtransz0(2,1:2),'all');
mtransz0(1,3:4) = mtransz0(1,3:4)/sum(mtransz0(1,3:4),'all');
mtransz0(2,3:4) = mtransz0(2,3:4)/sum(mtransz0(2,3:4),'all');
mtransz0(3,1:2) = mtransz0(3,1:2)/sum(mtransz0(3,1:2),'all');
mtransz0(4,1:2) = mtransz0(4,1:2)/sum(mtransz0(4,1:2),'all');
mtransz0(3,3:4) = mtransz0(3,3:4)/sum(mtransz0(3,3:4),'all');
mtransz0(4,3:4) = mtransz0(4,3:4)/sum(mtransz0(4,3:4),'all');

% krusell and smith (1998) assumes the idiosyncratic shock transition
% probability depends upon the aggregate shock realization. Thus, the 
% stationary equilibrium calculation needs an extra step to compute the
% stationary transition probability as follows:

%uu
mtransz(1,1)= mtransz0(1,1)*mtransA(1,1)...
             +mtransz0(1,3)*mtransA(1,2)...
             +mtransz0(3,1)*mtransA(2,1)...
             +mtransz0(3,3)*mtransA(2,2);
%ue
mtransz(1,2)= mtransz0(1,2)*mtransA(1,1)...
             +mtransz0(1,4)*mtransA(1,2)...
             +mtransz0(3,2)*mtransA(2,1)...
             +mtransz0(3,4)*mtransA(2,2);
%eu
mtransz(2,1)= mtransz0(2,1)*mtransA(1,1)...
             +mtransz0(2,3)*mtransA(1,2)...
             +mtransz0(4,1)*mtransA(2,1)...
             +mtransz0(4,3)*mtransA(2,2);
%ee
mtransz(2,2)= mtransz0(2,2)*mtransA(1,1)...
             +mtransz0(2,4)*mtransA(1,2)...
             +mtransz0(4,2)*mtransA(2,1)...
             +mtransz0(4,4)*mtransA(2,2);


% Wealth grid
pnumgrida   = 100;

% Finer grid near smaller wealth (Maliar, Maliar, and Valli, 2010)
vgridamin   = 0;
vgridamax   = 300;
x           = linspace(0,0.5,pnumgrida);
y           = x.^5/max(x.^5);
vgrida     = vgridamin+(vgridamax-vgridamin)*y;


%============
% Main Code %
%============

% 1. Compute stationary z distribution

% Eigenvector method
[vL,~] = eigs(mtransz',1);
vL = vL/sum(vL);
disp(vL);

% 2. Calculate aggregate supply of labour

supplyL = vgridz*vL;

% 3. Guess aggregate K -> prices {r,w}

% Initial guess
% Setup BEFORE the GE while loop:
K_ss = 30;
K_lo = K_ss * 0.7;
K_hi = K_ss * 1.3;
K    = 0.5*(K_lo + K_hi);
currentdist = ones(pnumgrida,pnumgridz)/(pnumgrida*pnumgridz);

% This is the outer loop (GE)
error2 = 10;
tol_ge = 1e-8;
pnumiter_ge = 1;

% Initialise
V = zeros(pnumgrida, pnumgridz);
pol_aprime = zeros(pnumgrida, pnumgridz);

while error2 > tol_ge

% Given K, all the prices are known.
r   = pAlpha*(K/supplyL)^(pAlpha-1)-pDelta;
mmu = r+pDelta;
w   = (1-pAlpha)*(K/supplyL)^(pAlpha);

% 4. VFI

% Pre-compute income (cash-on-hand)
gY = zeros(pnumgrida, pnumgridz);
for iz = 1:pnumgridz
    gY(:, iz) = (1+r)*vgrida(:) + w*vgridz(iz);
end

errorvfi = 10;
tolvfi   = 1e-7;
iter     = 0;

while errorvfi > tolvfi && iter < 2000
    iter   = iter + 1;
    V_prev = V;
    
    % Continuation value (discounted, expected over z')
    CV = pBeta * V * mtransz';   % pnumgrida x pnumgridz
    
    for iz = 1:pnumgridz
        CV_col = CV(:, iz);
        
        for ia = 1:pnumgrida
            coh = (1+r)*vgrida(ia) + w*vgridz(iz);
            
            % Step 1: grid search
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
            
            % Step 2: refine off-grid (this uses INTERPOLATION inside obj)
            ilo  = max(ap_best - 1, 1);
            ihi  = min(ap_best + 1, pnumgrida);
            a_lo = vgrida(ilo);
            a_hi = min(vgrida(ihi), coh - 1e-8);
            
            if a_hi > a_lo + 1e-10
                % interp1 inside obj estimates CV at any a' in [a_lo, a_hi]
                obj = @(ap) -(log(coh - ap) + interp1(vgrida, CV_col, ap, 'linear'));
                [aprime_opt, neg_val] = fminbnd(obj, a_lo, a_hi);
                pol_aprime(ia, iz) = aprime_opt;
                V(ia, iz)          = -neg_val;
            else
                pol_aprime(ia, iz) = vgrida(ap_best);
                V(ia, iz)          = V_best;
            end
        end
    end
    
    errorvfi = max(abs(V(:) - V_prev(:)));
end

fprintf('  VFI done in %d iters | pol_aprime range [%.3f, %.3f]\n', ...
    iter, min(pol_aprime(:)), max(pol_aprime(:)));

% 5. Compute the invariant distribution, Phi(a,z).
% NOTE: Make sure "currentdist = ones(pnumgrida,pnumgridz)/(pnumgrida*pnumgridz);" 
% is moved OUTSIDE and BEFORE the main GE while-loop for warm-starting!

errorhist = 10;
tolhist = 1e-8;
currentdist = ones(pnumgrida,pnumgridz)/(pnumgrida*pnumgridz);


while errorhist > tolhist
    nextdist = zeros(pnumgrida, pnumgridz);
    for iz = 1:pnumgridz
        for ia = 1:pnumgrida
            mass = currentdist(ia, iz);
            if mass < 1e-12; continue; end % Optimization: skip empty cells
            
            nexta = pol_aprime(ia, iz);
            nexta = min(max(nexta, vgrida(1)), vgrida(end));
            
            lb = sum(vgrida < nexta); 
            if lb < 1, lb = 1; end
            if lb >= pnumgrida, lb = pnumgrida - 1; end
            ub = lb + 1;
            
            weightlb = (vgrida(ub) - nexta) / (vgrida(ub) - vgrida(lb));
            weightlb = min(max(weightlb, 0), 1);
            weightub = 1 - weightlb;
            
            % Correct allocation across future idiosyncratic states
            for futureiz = 1:pnumgridz
                nextdist(lb, futureiz) = nextdist(lb, futureiz) + mass * mtransz(iz, futureiz) * weightlb;
                nextdist(ub, futureiz) = nextdist(ub, futureiz) + mass * mtransz(iz, futureiz) * weightub;
            end
        end
    end
    
    errorhist = max(abs(nextdist - currentdist), [], "all");
    currentdist = nextdist;
end
fprintf('  Invariant Distribution Found!\n');

% 6. Compute equilibrium allocations. We will close the GE loop.

% Recall that impK' = int{a'(a,z,S)}dPhi
%             K = int{a}dPhi(a,z)

marginaldista = sum(currentdist, 2); % Recall Phi(a,z), this sums rowwise. Gets the mass for each row (a)
endoK = vgrida(:)' * marginaldista(:);
%endoK = vgrida * marginaldista; % This is our aggregate K!

% 7. Check error2 and update

error2 = abs(endoK - K);
if endoK > K
    K_lo = K;
else
    K_hi = K;
end
K = 0.5*(K_lo + K_hi);

fprintf('GE %d | K=%.4f endoK=%.4f bracket=[%.4f, %.4f] err=%.2e\n', ...
    pnumiter_ge, K, endoK, K_lo, K_hi, error2);

pnumiter_ge = pnumiter_ge+1;

end

%================
% Store Results %
%================

% Equilibrium prices and aggregates
results.K        = K;
results.r        = r;
results.w        = w;
results.supplyL  = supplyL;
results.K_ss     = K_ss;

% Value function and policies
results.V          = V;
results.pol_aprime = pol_aprime;
results.pol_c      = (1+r)*repmat(vgrida(:), 1, pnumgridz) + w*repmat(vgridz, pnumgrida, 1) - pol_aprime;

% Stationary distribution
results.dist          = currentdist;
results.marginaldista = sum(currentdist, 2);
results.marginaldistz = sum(currentdist, 1);

% Save to disk for use in aggregate-uncertainty solver
save('ks1998_ss.mat', 'results');
fprintf('\nResults saved to ks1998_ss.mat\n');
fprintf('K* = %.4f | r* = %.5f | w* = %.4f\n', K, r, w);