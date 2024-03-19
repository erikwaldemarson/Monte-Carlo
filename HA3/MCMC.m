function [lambda,theta, t, accept_rate] = MCMC(tau, psi, rho, d, burn_in, M)
    %MCMC algorithm
    N = M + burn_in;
    N_accept = 0;
    
    t_0 = 1658;
    t_f = 1980;
    
    %initial values
    t_initial = linspace(t_0, t_f, d + 1)'; %end / breakpoints
    theta_initial = gamrnd(2, 1/psi);
    lambda_initial = gamrnd(2, 1/theta_initial, d, 1);

    %markov chains
    t = zeros(d+1,N); %chain of d+1 sized vectors
    lambda = zeros(d,N); %chains of d sized vector
    theta = zeros(1,N); %chain of thetas

    %start of MCMC 
    t(:,1) = t_initial; 
    lambda(:,1) = lambda_initial;
    theta(:,1) = theta_initial;

    idx = 1;

    for k = 1:(N-1)
        idx = idx + 1;
    
        if idx == d + 1
            idx = 2;
        end
        
        R = rho*(t(idx + 1, k) - t(idx - 1, k));
    
        %Random walk proposal: All at once
        epsilon = R*(2*rand(1) - 1); %Random walk U(-R, R)
        
        %proposal, r(z | t_k) = symmetric
        t_proposal = t(:,k); 
        t_proposal(idx) = t(idx, k) + epsilon; 

        %changes to order of breakpoints should have zero probability
        if ~issorted(t_proposal)
            alpha = 0; 
        else
            ratio = ln_f(t_proposal,lambda(:,k), tau) - ln_f(t(:,k), lambda(:,k), tau);
            alpha = min(1, exp(ratio));
        end

        u = rand();
        if u <= alpha
            t(:,k+1) = t_proposal;
            N_accept = N_accept + 1;
        else
            t(:,k+1)  = t(:,k);
        end

        %Gibbs step #1: theta_k+1
        a = 2*d + 2;
        b = psi + sum(lambda(:,k));
        theta(k+1) = gamrnd(a, 1/b);

        %Gibbs step #2: lambda_k+1
        for i = 1:d
            n_disasters = sum(tau >= t(i,k+1) & tau < t(i+1,k+1)); %number of distasters in interval

            a = n_disasters + 2;
            b = t(i+1,k+1) - t(i, k+1) + theta(k+1);
            lambda(i,k+1) = gamrnd(a, 1/b);
        end

    end

    accept_rate = N_accept/N;
    %end of MCMC
end

