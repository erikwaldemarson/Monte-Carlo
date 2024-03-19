%% load
load coal_mine_disasters.mat

t_0 = 1658;
t_f = 1980;

tau; %dates of disaster
indices = 1:numel(tau);

%plot(tau, indices)

%% Values

%initial values
burn_in = 500; %chain length before reaching stationary dist.
M = 5000; %samples in markov chain used for estimation
N = M + burn_in; %number of samples for each block


psi = 0.8;

d = 5; %number of intensities = nr. breakpoints + 1
% d > 1

%% 1.2 Hybrid MCMC

rho = 10;

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

N_accept = 0;
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

%end of MCMC
accept_rate = N_accept/N;

%% plots
figure

for i = 1:d
    subplot(d, 1, i); % Create subplot for each row
    plot(lambda(i,:));
    xlabel('n');
    ylabel('\lambda');

    if i == 1
        title('Plots of intensities \lambda for each block');
    end

    subtitle(['Block ', num2str(i)]);
    grid on;
end

figure   
plot(theta);
xlabel('n');
ylabel('\theta');
title(['Plot of \theta for d = ', num2str(d)]);
    
%%Using stationary from here on

t = t(:,burn_in+1:N);
lambda = lambda(:,burn_in+1:N);
theta = theta(burn_in+1:N);

%%Histograms 
for i = 1:d
    figure
    histogram(lambda(i, :));
    xlabel('\lambda');
    ylabel('Frequency');
    title(['Histogram of \lambda for Block ', num2str(i)]);
    grid on;
end

figure
histogram(theta)
xlabel('\theta');
ylabel('Frequency');
title('Histogram of \theta');


%% 1.3 Max-likelihood approach
burn_in = 500; %chain length before reaching stationary dist.
M = 5000; %samples in markov chain used for estimation
N = M + burn_in; %number of samples for each block
psi = 0.8;

rho_vec = [0.005, 0.01, 0.01, 0.03];

accept_rates = zeros(1,4);
likelihoods = zeros(1,4);

for i = 1:4
    d = i + 1; %2, 3, 4, 5
    rho = rho_vec(i);

    [lambda,theta, t, accept_rate] = MCMC(tau, psi, rho, d, burn_in, M);
    accept_rates(i) = accept_rate;

    test_cdf = [theta', gamcdf(theta,2,1/psi)'];


    mean_lambda = mean(lambda, 2);
    mean_t = mean(t, 2);

    likelihoods(i) = ln_f2(tau, mean_lambda, mean_t);
end




%% Estimates 

mean_lambda = mean(lambda, 2);
mean_theta = mean(theta);
mean_t = mean(t, 2);

% alpha = 0.01;
% z_alpha = norminv(1-alpha);
% z_alpha2 = norminv(1-alpha/2);

% L_lambda = mean_lambda - z_alpha.*std_lambda/sqrt(M);
% I_theta = [mean_theta - z_alpha.*std_theta/sqrt(M), mean_theta + z_alpha.*std_theta/sqrt(M)];

std_lambda = sqrt(var(lambda, 0, 2));
std_theta = sqrt(var(theta));

%% 1.4 Effect of prior
burn_in = 500; %chain length before reaching stationary dist.
M = 5000; %samples in markov chain used for estimation
N = M + burn_in; %number of samples for each block

d = 2;

rho = 0.005;

n = 10;
psi_vec = linspace(0.1, 5, n);
mean_theta = zeros(1, n);
mean_lambda = zeros(d,n);

for i = 1:n
    [lambda,theta, t, accept_rate] = MCMC(tau, psi_vec(i), rho, d, burn_in, M);
    mean_theta(i) = mean(theta);
    mean_lambda(:, i) = mean(lambda, 2);
    accept_rate

end

figure 
plot(psi_vec, mean_theta)
xlabel('\Psi')
ylabel('mean(\theta)')
title('mean(\theta) for different \Psi')

for i = 1:d
    figure
    plot(psi_vec, mean_lambda(i, :))
    xlabel('\Psi')
    ylabel('mean(\lambda)')
    title(['mean(\lambda) for Block ', num2str(i), ' for different \Psi'])
end





%% 1.5 Effect of rho on Mixing for theta
burn_in = 500; %chain length before reaching stationary dist.
M = 5000; %samples in markov chain used for estimation
N = M + burn_in; %number of samples for each block
psi = 0.8;
d = 2;

n = 10;
rho_vec = linspace(0.005, 10, n);

figure
for i = 1:n
    [lambda, theta, t, accept_rate] = MCMC(tau, psi, rho_vec(i), d, burn_in, M);

    [acf, lags] = autocorr(theta, NumLags=100);
   
    plot(lags, acf)
    hold on
    legend_strings{i} = sprintf('rho = %d', rho_vec(i));

end
xlabel('l');
ylabel('r(l) / \theta');
title('Plot of correlation function for \theta');
legend(legend_strings, 'Location', 'best');

hold off

%% 1.5 mixing for lambda

burn_in = 500; %chain length before reaching stationary dist.
M = 5000; %samples in markov chain used for estimation
N = M + burn_in; %number of samples for each block
psi = 0.8;
d = 2;

n = 10;
rho_vec = linspace(0.005, 10, n);



legend_strings = cell(1, n);
for i = 1:n
    [lambda, theta, t, accept_rate] = MCMC(tau, psi, rho_vec(i), d, burn_in, M);

    [acf, lags] = autocorr(lambda(1,:), NumLags=100); %assuming stationary process

    plot(lags, acf);
    xlabel('l');
    ylabel('r(l) / \lambda');
    legend_strings{i} = sprintf('rho = %d', rho_vec(i));
    
    hold on    
    
end
title('Plot of correlation function for \lambda_1');
legend(legend_strings, 'Location', 'best');
hold off




%% 
figure

for i = 1:n
    [lambda, theta, t, accept_rate] = MCMC(tau, psi, rho_vec(i), d, burn_in, M);

    for i = 1:d
        subplot(d, 1, i); % Create subplot for each row
        [acf, lags] = autocorr(lambda(i,:), NumLags=100); %assuming stationary process
    
        plot(lags, acf);
        xlabel('l');
        ylabel('r(l) / \lambda');
        subtitle(['Block ', num2str(i)]);
        legend_strings{i} = sprintf('Value: %d', rho_vec(i));

        if i == 1
            title('Plots of correlation function for \lambda for each block');
        end
        
        grid on;
        hold on
    end    
    
end
legend(legend_strings, 'Location', 'best');
hold off

%% 2.1
load atlantic.txt

f_inv = @(u, mu, beta) mu - beta.*log(-log(u));

%% 2.2
n = length(atlantic);
B = 200;
[beta, mu] = est_gumbel(atlantic);

boot_mu = zeros(1,B);
boot_beta = zeros(1,B);

for b = 1:B % bootstrap
    u = rand(1,n);
    y_boot = f_inv(u, mu, beta);
    [boot_beta(b), boot_mu(b)] = est_gumbel(y_boot);
end
delta_beta = sort(boot_beta - beta); % sorting to obtain quantiles
delta_mu = sort(boot_mu - mu); % sorting to obtain quantiles

alpha = 0.01; % CB level

L_beta = beta - delta_beta(ceil((1 - alpha/2)*B)); % forming CB beta
U_beta = beta - delta_beta(ceil(alpha*B/2));

L_mu = mu - delta_mu(ceil((1 - alpha/2)*B)); % forming CB mu
U_mu = mu - delta_mu(ceil(alpha*B/2));

I_beta = [L_beta, U_beta]
I_mu = [L_mu, U_mu]

%% 2.3

T = 3*14*100;
n = length(atlantic);
B = 200;

[beta, mu] = est_gumbel(atlantic);
tau_hat = f_inv(1-1/T, mu, beta);

boot = zeros(1,B);

for b = 1:B % bootstrap
    u = rand(1,n);
    y_boot = f_inv(u, mu, beta);
    [beta_boot, mu_boot] = est_gumbel(y_boot);
    
    boot(b) = f_inv(1-1/T, mu_boot, beta_boot);
end

delta = sort(boot - tau_hat); % sorting to obtain quantiles
alpha = 0.01; % CB level

U = tau_hat - delta(ceil(alpha*B))


