clc;
clear all;

rng(1)                       % set seed for random generator 
load('powercurve_D240.mat'); % load the power curve function


%% Initiate parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = 0:0.01:35;                                                      % wind speeds from 0 to 30 m/s with grid of 0.01 steps

month = [ "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December" ];
lambda = [ 11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10.0 10.9 11.7 11.7 ]; % scale parameter
k = [ 2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0 ];            % shape parameter


%% LOAD POWER CURVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = P(v);                    % get generated power

figure
plot(v, p./1000)             % scaling to get kW
ylim([0,16000])
title('Power curve for 2020ATB NREL Reference 15MW 240 (15 MW)')
xlabel('Wind speed v [m/s]')
ylabel('Output power P(v) [kW]')


%% PLOT WIND SPEED DISTRIBUTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 
for i = 1:12
    plot( v, wblpdf( v, lambda(i), k(i) ), 'DisplayName', sprintf('lamda=%G, k=%G, %s', lambda(i), k(i), month(i) ) )
    hold on
end

title('Wind speed distribution')
ylabel('pdf (Weibull)')
xlabel('Wind speed v [m/s]')
legend('show')


%% Question 2.a.1: STANDARD MONTE CARLO SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw N random samples from Weibull distribution and compute sample mean and standard deviation 

N = 1000;
z_score = norminv(0.995, 0, 1);   

tau = zeros(1,12);
stdev = zeros(1,12);
% I = zeros(2, 12);
for i = 1:12  
    windSample = wblrnd( lambda(i), k(i), [N 1] ); % rand generator of weibull 
    powerSample = P(windSample);
    tau(i) = mean(powerSample);
    stdev(i) = std(powerSample);
    I = [tau(i) - z_score * stdev(i)./sqrt(N), tau(i) + z_score * stdev(i)./sqrt(N)]/10^6
    % I(1, i) = (tau(i) + z_score * stdev(i)./sqrt(N) )./10^6;
    % I(2, i) = (tau(i) - z_score * stdev(i)./sqrt(N) )./10^6;
end

figure
plot( tau./10^6, 'b--' )                              % scaling to get MW
hold on 
plot( ( tau + z_score * stdev./sqrt(N) )./10^6, 'r' ) % scaling to get MW
hold on 
plot( ( tau - z_score * stdev./sqrt(N) )./10^6, 'r' ) % scaling to get MW
title( 'Standard Monte Carlo, N=1000, 99% confidence interval' )
xlabel( 'Month' )
ylabel( 'E[P(v)] / MW' )
legend( 'mean', 'bounds' )
xlim( [1 12] )

tau_std_standardMC = stdev./sqrt(N)./1000;                         % std of tau in kW for N simulations


%% Question 2.a.2: CONDITIONAL INVERSE MONTE CARLO SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw N random samples from Uniform distribution and plug into conditional 
% inverse cumjulative distribution function
%

N = 1000;
z_score = norminv(0.995, 0, 1);   
a = 4;
b = 25;

tau = zeros(1,12);
stdev = zeros(1,12);
eta = zeros(1,12);

for i = 1:12
    u = rand(1,N);
    Fa = wblcdf( a, lambda(i), k(i) );                           % cdf
    Fb = wblcdf( b, lambda(i), k(i) );                           % cdf
    eta(i) = Fb - Fa;
    windSample = wblinv( u .* (Fb - Fa) + Fa, lambda(i), k(i) ); % inv cdf
    powerSample = P(windSample);
    tau(i) = mean(powerSample);
    stdev(i) = std(powerSample);
    eta(i);
    I = eta(i)*[tau(i) - z_score * stdev(i)./sqrt(N), tau(i) + z_score * stdev(i)./sqrt(N)]/10^6
end

figure
plot( tau./10^6, 'b--' )                              % scaling to get MW
hold on 
plot( ( tau + z_score * stdev./sqrt(N) )./10^6, 'r' ) % scaling to get MW
hold on 
plot( ( tau - z_score * stdev./sqrt(N) )./10^6, 'r' ) % scaling to get MW
title( 'Truncated Monte Carlo, N=1000, 99% confidence interval' )
xlabel( 'Month' )
ylabel( 'E[P(v)] / MW' )
legend( 'mean', 'bounds' )
xlim( [1 12] )


%% Question 2.b: CONTROLE VARIATE MONTE CARLO SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw N random samples from Uniform distribution and plug into conditional 
% inverse cumjulative distribution function
%
% P_{tot}(V) = 0.5*rho*pi*(d^2)*(v^3)/4 -> Total power produced by wind mill [W]
% P(V) = alpha*P_{tot}(V)               -> Power produced in practice scaled down by the power coefficient alpha
%
% gamma = 0.5*rho*pi*(d^2)/4
% alpha = UNKNOWN!                      -> must be estimated or calculated explicitly, we cannot calculate explicitly since we do not know alpha!
% beta =  - Cov(V, P(V)) / Var(V)       -> must be estimated or calculated explicitly
%
% CV = P(V) + beta * ( V - E[V] )       -> Controle variate which we want to estimate
%

N = 1000;
z_score = norminv(0.995, 0, 1);                    % two sided 99% confidence interval


tau = zeros(1,12);
stdev = zeros(1,12);
beta = zeros(1,12);

for i = 1:12   
    windSample = wblrnd( lambda(i), k(i), [N 1] ); % rand generator from Weibull 
    powerSample = P(windSample);                   % P(v)
   
    %windSample = windSample.^3;                                   % transfromation of wind to get covariance different from zero 
    %var_V = E_V(lambda(i), k(i), 6) - E_V(lambda(i), k(i), 3).^2; %explicit variance calculation V^3
    
    var_V = E_V(lambda(i), k(i), 2) - E_V(lambda(i), k(i), 1).^2; % explicit variance calculation V 
    %var_V = var(windSample);                                      % sample estimate -> Biased beta!!

    c = cov(powerSample, windSample);
    beta(i) = - c(2,1) / var_V;

    %controlvariate = powerSample + beta(i) * (windSample - E_V(lambda(i), k(i), 3)); % variate = V^3
    controlvariate = powerSample + beta(i) * (windSample - E_V(lambda(i), k(i), 1)); % variate = V
    

    tau(i) = mean(controlvariate);
    stdev(i) = std(controlvariate);

    I = [tau(i) - z_score * stdev(i)./sqrt(N), tau(i) + z_score * stdev(i)./sqrt(N)]/10^6
end

figure
plot( tau./10^6, 'b--' )                              % scaling to get MW
hold on 
plot( ( tau + z_score * stdev./sqrt(N) )./10^6, 'r' ) % scaling to get MW
hold on 
plot( ( tau - z_score * stdev./sqrt(N) )./10^6, 'r' ) % scaling to get MW
title( 'Control Variate Monte Carlo, N=1000, 99% confidence interval' )
xlabel( 'Month' )
ylabel( 'E[P(v)] / MW' )
legend( 'mean', 'bounds' )
xlim( [1 12] )




%% Question 2.c.1: IMPORTANCE SAMPLING MONTE CARLO SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visual inspection of optimal choice of g(x)
% P(v)f(v) = 0 and g(v) = 0 must be true for the same value of v! -> Condition for choice of g(v)
% Also the distribution g(v) must follow the distribution of P(v)f(v) fairly well!
% Also the distribution g(v) must be easy to draw random samples from!
% Also the ration P(v)f(v) / g(v) must be bounded for all v! 
%
% -> We test with the Rayleigh distribution since it is realetd to the Weibull distribution!
% -> Optimal value of sigma?
%

% Sigma must be > 0 but when sigma come close to 0 then the ratio goes to infinity!
% Good choice of sigma seems to be 6, 8, 10, 12!
month_i = 3;
sigma = 8;             % adjust to find optimal value
alpha = 6.2;
beta = 2;
scaleFacProd = 12000000; % adjust to find optimal value
scaleFacPow = 250000000; % adjust to find optimal value

powerSample = P(v);
f = wblpdf(v, lambda(month_i), k(month_i));
%g = raylpdf( v, sigma );
g = gampdf(v, alpha, beta);
prod = powerSample .* f';

figure 
subplot(2,1,1)
plot( v, g, 'DisplayName', sprintf('g(v), sigma=%G', sigma ) )
hold on
plot( v, f, 'DisplayName', sprintf('f(v), lambda=%G, k=%G', lambda(month_i), k(month_i) ) )
hold on
plot( v, powerSample ./ scaleFacPow, 'DisplayName', sprintf('P(v) / %G', scaleFacPow) )
hold on
plot(v, prod ./ scaleFacProd, 'DisplayName', sprintf('P(v)f(v) / %G', scaleFacProd) )
title( sprintf('g(x)=Rayleigh - %s', month(month_i)) )
ylabel('Distribution')
xlabel('Wind speed v [m/s]')
legend('show')

subplot(2,1,2)
plot(v, (prod ./ scaleFacProd) ./g')
title( 'Ratio: P(v)f(v) / g(v)' )
xlabel('Wind speed v [m/s]')


%% Question 2.c.2: INSTRUMENTAL DISTRIBUTION MONTE CARLO SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 1000;
z_score = norminv(0.995, 0, 1);   
alpha = 6.2;
beta = 2;

tau = zeros(1,12);
stdev = zeros(1,12);
instrumentalSample = zeros(N,12);

for i = 1:12  
    windSample = gamrnd(alpha, beta, [N 1] ); % rand generator of gamma
    g = gampdf(windSample, alpha, beta);
    f = wblpdf(windSample, lambda(i), k(i));
    powerSample = P(windSample);
    
    %instrumentalSample(:, i) = powerSample .* (f./g) ./ (14*max(powerSample .* f));
    instrumentalSample(:, i) = powerSample .* (f./g);
    
    %tau(i) = mean(instrumentalSample(:, i)) .* (14*max(powerSample .* f));
    %stdev(i) = std(instrumentalSample(:, i)) .* (14*max(powerSample .* f));

    tau(i) = mean(instrumentalSample(:, i)) ;
    stdev(i) = std(instrumentalSample(:, i)) ;

     I = [tau(i) - z_score * stdev(i)./sqrt(N), tau(i) + z_score * stdev(i)./sqrt(N)]/10^6
end

figure
plot( tau./10^6, 'b--' )                              % scaling to get MW
hold on 
plot( ( tau + z_score * stdev./sqrt(N) )./10^6, 'r' ) % scaling to get MW
hold on 
plot( ( tau - z_score * stdev./sqrt(N) )./10^6, 'r' ) % scaling to get MW
title( 'Importance sampling Monte Carlo, N=1000, 99% confidence interval' )
xlabel( 'Month' )
ylabel( 'E[P(v)] / MW' )
legend( 'mean', 'bounds' )
xlim( [1 12] )




%% Question 2.d: ANTITHETIC VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 1000;
z_score = norminv(0.995, 0, 1);   


tau = zeros(1,12);
stdev = zeros(1,12);
correl = zeros(1,12);

for i = 1:12
    u = rand(1,N);
    windSample = wblinv( u, lambda(i), k(i) ); % inv cdf
    windSampleStar = wblinv( ones(1,N) - u, lambda(i), k(i) ); % inv cdf
    powerSample = P(windSample);
    powerSampleStar = P(windSampleStar);
    tau(i) = mean( (powerSample + powerSampleStar) ./ 2);
    stdev(i) = std( (powerSample + powerSampleStar) ./ 2);
    correl(i) = corr(powerSampleStar, powerSample);
    I = [tau(i) - z_score * stdev(i)./sqrt(N), tau(i) + z_score * stdev(i)./sqrt(N)]/10^6;
end

figure
plot( tau./1000, 'b--' )                              % scaling to get kW
hold on 
plot( ( tau + z_score * stdev./sqrt(N) )./1000, 'r' ) % scaling to get kW
hold on 
plot( ( tau - z_score * stdev./sqrt(N) )./1000, 'r' ) % scaling to get kW
title( 'Antithetic Variable Monte Carlo, N=1000, 99% confidence interval' )
xlabel( 'Month' )
ylabel( 'E[P(v)] [kW]' )
legend( 'mean' )
xlim( [1 12] )

tau_std_antithetic = stdev./sqrt(N)./1000;                         % std of tau in kW for N simulations


%% Question 2.e: Probability of positive power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prob = zeros(1,12);

for i = 1:12
    prob(i) = wblcdf( 25, lambda(i), k(i) ) - wblcdf( 4, lambda(i), k(i) ); % inv cdf
end

figure
plot( prob )                              % scaling to get kW


%% Question 2.f: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 1000;
z_score = norminv(0.995, 0, 1);   
theta = 0.5*1.225*pi*(240^2)/4;                                  % constant from P_tot function

tau = zeros(1,12);
stdev = zeros(1,12);
powerTot = zeros(1,12);
ratio = zeros(1,12);

for i = 1:12
    u = rand(1,N);
    windSample = wblinv( u, lambda(i), k(i) ); % inv cdf
    windSampleStar = wblinv( ones(1,N) - u, lambda(i), k(i) ); % inv cdf
    powerSample = P(windSample);
    powerSampleStar = P(windSampleStar);
    
    tau(i) = mean( (powerSample + powerSampleStar) ./ 2);
    powerTot(i) = theta * E_V( lambda(i), k(i), 3);

    stdev(i) = std( (powerSample + powerSampleStar) ./ 2) / powerTot(i);
    
    ratio(i) = tau(i)/powerTot(i);

    I = [ratio(i) - z_score * stdev(i)./sqrt(N), ratio(i) + z_score * stdev(i)./sqrt(N)]*100
end



figure
plot( ratio .* 100, 'b--' )                              % scaling to get kW
hold on 
plot( ( ratio + z_score * stdev./sqrt(N) ).*100, 'r' ) % scaling to get kW
hold on 
plot( ( ratio - z_score * stdev./sqrt(N)).*100, 'r' ) % scaling to get kW
title( 'Power coefficient with antithetic variable' )
xlabel( 'Month' )
ylabel( '%' )
legend( 'mean' )
xlim( [1 12] )

%ratio_std = stdev./sqrt(N).* 100                        % std of tau in kW for N simulations


%% Question 2.g: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From the results we see that the spot is a good place to plant the plant 

% capacity factor: Time with power production / total time
N = 1000;

tau = zeros(1,12);

for i = 1:12
    u = rand(1,N);
    windSample = wblinv( u, lambda(i), k(i) ); % inv cdf
    windSampleStar = wblinv( ones(1,N) - u, lambda(i), k(i) ); % inv cdf
    powerSample = P(windSample);
    powerSampleStar = P(windSampleStar);
    tau(i) = mean( (powerSample + powerSampleStar) ./ 2);
end

sum(tau) ./ (15000000*12)

% availability factor: Time with power production / total time
prob = zeros(1,12);

for i = 1:12
    prob(i) = wblcdf( 25, lambda(i), k(i) ) - wblcdf( 4, lambda(i), k(i) ); % inv cdf
end

mean(prob)


%% Question 3.a: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1.95;
lambda = 10.05;

N = 1000;


sigma = 12;
alpha = 6.2;
beta = 2;


u = rand(1,N);
windSample = wblinv( u, lambda, k); % inv cdf
windSampleStar = wblinv( ones(1,N) - u, lambda, k ); % inv cdf
powerSample = P(windSample);
powerSampleStar = P(windSampleStar);
tau = mean( (powerSample + powerSampleStar) ./ 2)
stdev = std( (powerSample + powerSampleStar) ./ 2)


tau_sum = 2*tau;
stdev_sum = 2*stdev;

%% Functions: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1.95;
lambda = 10.05;
a = 0.638;
p = 3;
q = 1.5;

 
factor = @(x) (1-wblcdf(x, lambda, k).^p ).^(q-1).*(wblcdf(x, lambda, k).^(p).*(1+p*q)-1);
factor2 = @(x) (1-wblcdf(x, lambda, k).^p ).^(q);

wb2Dpdf = @(x, y) wblpdf(x, lambda, k).*wblpdf(y, lambda, k).*(1+a*factor(x).*factor(y));
wb2Dcdf = @(x, y) wblcdf(x, lambda, k).*wblcdf(y, lambda, k)*(1+a*factor2(x)*factor2(y));



%% 3b

n = 100;
gridV = linspace(0, 35, n);
[X,Y] = meshgrid(gridV',gridV');

mu = [11.5 11.5];
Sigma = [30 7;
        7 30];
f = mvnpdf([X(:), Y(:)], mu, Sigma);
    
% C = [20 7;
%     7 20];
% df = 3;
% 
% f = mvtpdf([X(:), Y(:)], C, df); 


prod = wb2Dpdf(X(:),Y(:)).*P(X(:)').*P(Y(:)') / (12.5^11);
       
%mean = sum / (length(gridV)  + length(gridV))

figure(1)
surf(X,Y,reshape(prod./f, n, n));
figure(2)
surf(X,Y,reshape(prod, n, n));
figure(3)
surf(X,Y,reshape(f, n, n));


%% Question 3.b: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau = 7.6*10^6; %3a

N = 1000;

mu = [11.5 11.5];
Sigma = [30 7;
        7 30];


windSample = mvnrnd(mu, Sigma, N);

f = wb2Dpdf(windSample(:,1),windSample(:, 2));
psi = P(windSample(:, 1)').* P(windSample(:, 2)');
g = mvnpdf(windSample, mu, Sigma);
instrumentalSample = f.*psi ./g;


tau_prod = mean(instrumentalSample)

CoV = tau_prod - tau^2


%% 3c

Var_sum = 2*stdev^2 + 2*CoV

%% 3d

N = 1000;
N_u = zeros(N, 1);
N_l = zeros(N,1);

z_score = norminv(0.995, 0, 1);                    % two sided 99% confidence interval

mu = [11.5 11.5];
Sigma = [30 7;
        7 30];
    
    
windSample = mvnrnd(mu, Sigma, N);

psi = zeros(N,1);

for i = 1:N
    prod = wb2Dpdf(windSample(i, 1),windSample(i, 2)) * 1;
    g = mvnpdf([windSample(i, 1) windSample(i, 2)] , mu, Sigma);
    
    if (P(windSample(i, 1)) + P(windSample(i, 2))) > 15*10^6
            N_u(i) = prod ./g;
    else
            N_l(i) = prod ./g; 
    end
end

tau_u = mean(N_u);
std_u = std(N_u);

tau_l = mean(N_l);
std_l = std(N_l);


I_u = [tau_u - z_score * std_u/sqrt(N), tau_u + z_score * std_u/sqrt(N)]
I_l = [tau_l - z_score * std_l/sqrt(N), tau_l + z_score * std_l/sqrt(N)]

total_prob = tau_u + tau_l


