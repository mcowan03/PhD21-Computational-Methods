% Solve for TCE %%
clear;
% Set the parameters
% Parameters
sigma = 1;
rho = 0.95;
sigma_epsilon = 0.009;
alpha = 0.33;
beta = 0.99;
delta = 0.025;
chi = 1;
mu = 0.60;


% Load the two steady states.
A_old = 1;
A_new = 1.1;
[SS1, eta] = ssSolver(A_old);
[SS2] = ssSolverNewA(A_new);

% Load K and N values for both states
K1 = SS1.K;
N1 = SS1.N;
K2 = SS2.K;
N2 = SS2.N;

% Load prices
r_ss2 = SS2.r;
w_ss2 = SS2.w;

% Set a time path
T = 200;

% Initialise
vA = 1.100 * ones(T,1);
vK = K2 * ones(T+1,1); % K path guess
vK(1) = K1;
vN = N2 * ones(T,1); % N path guess
vr = r_ss2 * ones(T,1); 
vw = w_ss2 * ones(T,1);
C_ss2 = r_ss2 * K2 + w_ss2 * N2;
vC = vr.*vK(1:T) + vw.*vN; % 
vCnew = vC; % C path guess
vNnew = vN; % We don't guess N path?
it = 0;
error1 = 1;

while error1 > 1e-7
    % Solving backwards - we are guessing on 
    vCprime = [vC(2:T); C_ss2]; % T length vector, corresponds to t+1
    vrprime = [vr(2:T); r_ss2]; % T length vector, corresponds to t+1
    vKprimeprime = [vK(3:T+1); K2]; % T length vector, corresponds to t+2
    vKprime = vK(2:T+1); % T length vector, corresponds to t+1
    vKthis = vK(1:T); % T length vector, corresponds to t
    vPsi1prime = mu/2 * (1 - (vKprimeprime./vKprime).^2); % Required for euler
    vPsi2 = mu * (1-vKprime./vKthis);
    tempC = (beta*(vCprime.^(-sigma).*(1+vrprime-vPsi1prime))./(1+vPsi2)).^(-1/sigma); % Euler. Updates for C.
    vNnew = (tempC.^(-sigma) .* vw / eta).^(chi); % Solves for N

    % Simulate Forwards
    vPsi = mu/2*(diff(vK)).^2./vK(1:T);
    vY = vA.*vK(1:T).^alpha.*vNnew.^(1-alpha);
    vrnew = alpha .* vA .* (vK(1:T)./vNnew).^(alpha-1) - delta;
    vwnew = (1-alpha) .* vA .* (vK(1:T)./vNnew).^(alpha);
    vInew = vY - tempC - vPsi;
    vKnew = vK;
    vKnew(2:(T+1)) = vInew + (1-delta) * vK(1:T);
    vCnew = vA.*vKnew(1:T).^alpha.*vNnew.^(1-alpha) - vInew + mu/2*(diff(vKnew)).^2./vKnew(1:T);

    error1 = max(abs([vK - vKnew; vC - vCnew; vw - vwnew; vr - vrnew]));

    if mod(it, 50) == 1
        % Live Plotting of K, N, C, and Y
        subplot(2, 2, 1);
        plot(1:(T+1), vK, 'b'); hold on;
        plot(1:(T+1), vKnew, 'r-.'); hold off;
        xlim([1,T+1]);
        title(['K (Iter: ', num2str(it), ')']); 
        legend("Guess","New");

        subplot(2, 2, 2);
        plot(1:T, vN, 'b'); hold on;
        plot(1:T, vNnew, 'r-.'); hold off;
        xlim([1,T]);
        title('N');
        legend("Guess","New");

        subplot(2, 2, 3);
        plot(1:T, vC, 'b'); hold on;
        plot(1:T, vCnew, 'r-.'); hold off;
        xlim([1,T]);
        title('C');
        legend("Guess","New");

        subplot(2, 2, 4);
        plot(1:T, vY, 'k-'); 
        xlim([1,T]);
        title('Y (Current Realized)');
        
        drawnow;
        
        Phrase = ['Iteration is in progress: ',num2str(it),'st iteration'];
        disp(Phrase);
        fprintf('Error: %.18f \n \n', error1);
    end



    % update
    vK = vKnew * 0.1 + vK * 0.9;
    vC = vCnew * 0.1 + vC * 0.9;
    vw = vwnew * 0.1 + vw * 0.9;
    vr = vrnew * 0.1 + vr * 0.9;

    it = it + 1;
end
