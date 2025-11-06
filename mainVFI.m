%% SOLVING AIYAGARI 1994 MODEL %%

% The algorithm:
% 1. Parameterize
% 2. Guess price, r
% 3. Solve HH problem - VFI
% 4. Compute stationary distribution
% 5. Update prices/aggr.s
% 6. Check convergence, if not repeat from 2.
verbose             = true; %reports on

% Parameterize

    beta = 0.96;
    delta = 0.08;
    alpha = 0.36;
    A = 1.0;

% Income process (z)
    m = 3.0; % m std. dev. either side of mu
    mu = 0;
    rho = 0.9;
    sigma = 0.2 * sqrt(1 - rho^2);
    Nz = 7; % number of grid points

% Wealth Grid
    a_lower = 0.0;
    a_upper = 150.0;
    Na = 200;
    Na2 = 400;

% Finer points at low wealth states: (Maliar, Maliar and Valli, 2010)
    x = linspace(0,0.5,Na);
    x2 = linspace(0,0.5,Na2);
    y = x.^5/max(x.^5);
    y2 = x2.^5/max(x2.^5);
    a_grid = a_lower+(a_upper-a_lower)*y;
    a_grid2 = a_lower+(a_upper-a_lower)*y2;

% Productivity

[logzgrid, P] = tauchen(rho, sigma, mu, Nz, m);
Pz = stationary_distribution(P);
zGrid = exp(logzgrid);


% Initialize

mVF = zeros(Na, Nz); 
mVFnew = mVF;
mCPol = zeros(Na, Nz);
mAPrimePol = zeros(Na, Nz);
%define the size of the distribution
currentDist = ones(Na2,Nz);
currentDist = currentDist/(Na2*Nz); %normalize




% Initial Guess

K = 6;
weightOld = 0.9;

[vL,~] = eigs(P',1);
vL = vL/sum(vL);
disp(vL);
supplyL = zGrid*vL

%% SOLVING THE MODEL %%
% Loop Settings
error2 = 10;
tol_ge = 1e-8;
numIterGE = 1;


% GE LOOP %
while error2 > tol_ge
    % guessing K
    r = alpha*(K/supplyL)^(alpha-1)-delta;
    w = (1-alpha)*(K/supplyL)^alpha;

% VFI %
error = 10;
numIterVFI = 1;

    while error > 1e-8
           for iz = 1:Nz
               z = zGrid(iz);
               % expected future value, given z
               expValue = mVF*P;
               expValue = expValue(iz,:);

               % We use monotonicity
               minWealth = a_lower;

               for ia = 1:Na

                   a = a_grid(ia);
                   budget = w*z+(1+r)*a;

                   % Optimal Saving
                   aprime = fnOptFAST(beta,budget,a_grid,expValue,Na,minWealth);

                   % Borrowing constraint
                   aprime(aprime<a_lower) = a_lower;

                   % Linear interpolation for off-the-grid aprime
                   a_lower = sum(a_grid<aprime);
                   a_lower(a_lower<=1) = 1;
                   a_lower(a_lower>=Na) = Na-1;
                   a_upper = a_lower+1;
                   weightLow = (a_grid(a_upper) - aprime)/(a_grid(a_upper)-a_grid(a_lower));
                   value = weightLow*expValue(a_lower)+(1-weightLow)*expValue(a_upper);

                   c = budget - aprime;

                   % Update
                   minWealth = aprime;
                   mVFnew(ia,iz) = log(c)+beta*value;
                   mCPol(ia,iz) = c;
                   mAPrimePol(ia,iz) = aprime;

               end
           end

       % The iteration
       error = max(max(max(abs(mVFnew-mVF))));
       mVF = mVFnew;

       
numIterVFI = numIterVFI+1;

    end %This ends the VFI loop

    %interpolation - wealth policy
    %=========================
    %define the size of interpolated policy function
    mAPrimePol2 = zeros(size(currentDist));
    %interpolate
    if Na == Na2
        mAPrimePol2 = mAPrimePol;
        else
        for iz = 1:Nz
        mAPrimePol2(:,iz) = interp1(a_grid,aprime(:,iz),a_grid2,"linear","extrap");
        end
    end

% Non Stochastic Simulation - Eigenvec method % Similar to KFE?
% Going to find a large transition probability matrix %

mPolicy = zeros(Na2*Nz, Na2*Nz);

vLocationCombineda = kron(1:Na2,ones(1,Nz));
vLocationCombinedz = kron(ones(size(a_grid2)),1:Nz);

for iLocation = 1:Nz*Na2

    ia = vLocationCombineda(iLocation);
    iz = vLocationCombinedz(iLocation);

    a = a_grid2(ia);
    nexta = mAPrimePol2(ia,iz);
    LB = sum(a_grid2<nexta);
    LB(LB<=0) = 1;
    LB(LB>=Na2) = Na2-1;
    UB = LB+1;
    weightLB = (a_grid2(UB) - nexta)/(a_grid2(UB)-a_grid2(LB));
    weightLB(weightLB<0) = 0;
    weightLB(weightLB>1) = 1;
    weightUB = 1-weightLB;

    for izprime = 1:Nz

        mPolicy(iLocation,:) = mPolicy(iLocation,:)+(vLocationCombineda==LB).*(vLocationCombinedz==izprime) * weightLB * P(iz,izprime);
        mPolicy(iLocation,:) = mPolicy(iLocation,:)+(vLocationCombineda==UB).*(vLocationCombinedz==izprime) * weightUB * P(iz,izprime);

    end

end

mPolicytrans = mPolicy';

[currentDist0,~] = eigs(mPolicytrans,1);%eig(mPolicy');
currentDist0 = currentDist0(:,1)/sum(currentDist0(:,1));
currentDist = zeros(Na2,Nz);

for iLocation = 1:Nz*Na2

    ia = vLocationCombineda(iLocation);
    iz = vLocationCombinedz(iLocation);

    currentDist(ia,iz) = currentDist0(vLocationCombineda==ia & vLocationCombinedz==iz);

end
currentDist(currentDist<0) = 0;

%endogenous aggregate allocation
marginalDista = sum(currentDist,2);
endoK = a_grid2*marginalDista;

%error and update
error2 = abs(endoK - K);
K = K.*weightOld+endoK.*(1-weightOld);

% report only spasmodically
if verbose == true && (floor((numIterGE-1)/50) == (numIterGe-1)/50) || error2<= tol_ge
%=========================  
% interim report
%=========================  

fprintf(' \n');
fprintf('market clearing results \n');
fprintf('max error: %.15f \n', error2);
fprintf('capital rent: %.15f \n', r);
fprintf('wage: %.15f \n', w);
fprintf('aggregate capital: %.15f \n', K);

% plot
close all;
figure;
plot(a_grid2,currentDist);
title("The wealth distributions for different labor endowments","fontsize",15)
saveas(gcf,'../figures/dist_ss.jpg');
hold off;

pause(0.01);
toc;

end

numIterGE = numIterGE + 1;

end