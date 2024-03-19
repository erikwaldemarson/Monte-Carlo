%HA2

%% 1.3 Naive SIS
d = 2;
N = 1000; %number of particles
n = 71;   %length of walk

%%Start of SIS
X = zeros(n, d); %initial X, each row represents a coordinate in the random walk
X(:,:, N) = 0; %each page represents a different particle
w = ones(N, n); %weights, particle per row and step per column

for k = 1:n
    for i = 1:N
        step = zeros(1, d); %random step
        idx = randi([1,d]); %choose random index
        step(1, idx) = 1;   %set to 1
        flip = randi([0,1]);%flip a coin
        
        if flip == 0        %if heads then multiply by minus -1
            step = -1*step;
        end
        
        X(k+1, :, i) = X(k, :, i) + step; %take a step
        
        w(i, k+1) = w(i, k)*4; %default
        
        if w(i,k) ~= 0
            for l = 1:k %seeing if we are self-avoiding
                if norm(X(k+1, :, i) - X(l,:, i)) == 0
                    w(i, k+1) = w(i, k)*0; %if not then set weight to 0
                    break
                end
            end
        end
     
    end
end
%%end of SIS
is
c_n = mean(w(:,n+1)); %c_n = average of weights


%% 1.4 Improved SIS
d = 2;
N = 100; %number of particles
n = 71;   %length of walk

%%Start of SIS
X = zeros(n, d); %initial X, each row represents a coordinate in the random walk
X(:,:, N) = 0; %each page represents a different particle
w = ones(N, n); %weights, particle per row and step per column

for k = 1:n
    for i = 1:N
        %self avoiding walk g

        isAllowed = ones(2,d); %0 if step is not allowed
        
        for j = 1:d
            for l = 1:k
                step = zeros(1,d);
                step(1,j) = 1;
                
                if X(l, :, i) == X(k, :, i) + step
                    isAllowed(1,j) = 0;
                end
                
                if X(l, :, i) == X(k, :, i) - step 
                    isAllowed(2,j) = 0;
                end
            end
        end
        
        nr_of_allowed = sum(isAllowed, 'all');
        set_of_allowed = cell(1, nr_of_allowed);
        
        %create set of all allowed steps
        x = 1;
        for l = 1:2
            for j = 1:length(isAllowed)
                if isAllowed(l,j) == 1
                    set_of_allowed{x} = [l,j];
                    x = x + 1;
                end
            end
        end
        
        if nr_of_allowed == 0 %if no steps are allowed then X_k+1 = X_k
            X(k+1, :, i) = X(k, :, i);
            w(i, k+1) = 0;
        else
            step = zeros(1, d); %random step
            pos = randi([1,nr_of_allowed]); %choose random index

            idx = set_of_allowed{pos}(2);
            
            step(idx) = 1;

            if set_of_allowed{pos}(1) == 2 
                step = -1*step;
            end
            
            X(k+1, :, i) = X(k, :, i) + step;
            w(i, k+1) = w(i, k)*nr_of_allowed; 

        end

     
    end
end
%%end of SIS

c_n = mean(w(:,n+1)) %c_n = average of weights


%% 1.5 SISR
d = 2;
N = 1000; %number of particles
n = 100;   %length of walk

%%Start of SISR
X = zeros(n, d); %initial X, each row represents a coordinate in the random walk
X(:,:, N) = 0; %each page represents a different particle
w = ones(N, n); %weights, particle per row and step per column

for k = 1:n
    for i = 1:N
        %selection
        z = 1:N;
        W = w(:,k)./sum(w(:,k));
        x_idx = randsample(z,1,true,W);
        
        %mutation
        X(:, :, i) = X(:, :, x_idx); 
        
        %self avoiding walk g
        isAllowed = ones(2,d); %0 if step is not allowed
        
        for j = 1:d
            for l = 1:k
                step = zeros(1,d);
                step(1,j) = 1;
                
                if X(l, :, i) == X(k, :, i) + step
                    isAllowed(1,j) = 0;
                end
                
                if X(l, :, i) == X(k, :, i) - step 
                    isAllowed(2,j) = 0;
                end
            end
        end
        
        nr_of_allowed = sum(isAllowed, 'all');
        set_of_allowed = cell(1, nr_of_allowed);
        
        %create set of all allowed steps
        x = 1;
        for l = 1:2
            for j = 1:length(isAllowed)
                if isAllowed(l,j) == 1
                    set_of_allowed{x} = [l,j];
                    x = x + 1;
                end
            end
        end
        
        if nr_of_allowed == 0 %if no steps are allowed then X_k+1 = X_k
            X(k+1, :, i) = X(k, :, i);
            w(i, k+1) = 0;
        else
            step = zeros(1, d); %random step
            pos = randi([1,nr_of_allowed]); %choose random index

            idx = set_of_allowed{pos}(2);
            
            step(idx) = 1;

            if set_of_allowed{pos}(1) == 2 
                step = -1*step;
            end
            
            X(k+1, :, i) = X(k, :, i) + step;
            w(i, k+1) = nr_of_allowed; 

        end
     
    end
end
%%end of SISR
c_n = 1;
for i = 2:n+1
    c_n = mean(w(:,i))*c_n; %c_n = average of weights
end

c_n

%% plotting SAW
plot(X(:,1,39),X(:,2,39))
xlabel('x')
ylabel('y')
title('Self-avoiding walk in Z^2 for n = 100')


%% 1.6 and 1.9
n_max = 50;
x = 1:n_max;
c_n = zeros(1,n_max);
%d = 2; %1.6
d = 5;  %1.9
N = 100;

for n = 1:n_max
    
    c_n(n) = SISR(d, n, N); 
    
end

%Regression model
y = log(c_n);

% Define your function
fun = @(coeffs, x) coeffs(1) + coeffs(2)*x + coeffs(3)*log(x);

% Initial guess for the coefficients
initial_guess = [0, 0, 0]; % [m, k, c]

% Perform nonlinear least squares fitting
coeffs_fit = lsqcurvefit(fun, initial_guess, x, y);

% Extract fitted coefficients
m_fit = coeffs_fit(1);
k_fit = coeffs_fit(2);
c_fit = coeffs_fit(3);

% Plot the original data and the fitted curve



figure;
plot(x, y, 'bo', 'DisplayName', 'Data');
hold on;
x_fit = linspace(min(x), max(x), 100);
y_fit = m_fit + k_fit*x_fit + c_fit*log(x_fit);
plot(x_fit, y_fit, 'r-', 'DisplayName', 'Fitted Curve');
xlabel('n');
%ylabel('ln c_n(2)'); %1.6
ylabel('ln c_n(5)'); %1.9
%title('Curve Fitting d = 2'); %1.6
title('Curve Fitting d = 5'); %1.9
legend('Location', 'best');
grid on;
hold off;

A_2 = exp(m_fit)
u_2 = exp(k_fit)
gamma_2 = 1 + c_fit

%% 1.9 continued
d = 5;

x = 2*d - 1 - 1/(2*d) - 3/(2*d)^2 - 16/(2*d)^3;

mu = [8.8160, 8.8264, 8.8231];

error = mu - x*ones(1,3);

error / (1/d^4)

%% 2.1
load population_2024.mat

%plot(X, Y,'*')

A = 0.8;
B = 3.8;
C = 0.6;
D = 0.99;
G = 0.8;
H = 1.25;

N = 10000;
n = 100;
tau = zeros(1,n+1); % vector of filter means
w = zeros(N,1);
p = @(x,y) unifpdf(y,G*x,H*x); % observation density, for weights
part = C + (D-C).*rand(N,1); % initialization
w = p(part,Y(1)); % weighting
tau(1) = sum(part.*w)/sum(w); % estimation

% conf int
[xx,I]=sort(part); % sort data
cw=cumsum(w(I))/sum(w); % cumulative normalized weightsum
% for sorted data
Ilower=find(cw>=0.025,1); % index for lower 2.5% quantile
Iupper=find(cw>=0.975,1); % index upper 2.5% quantile
taulower(1)=xx(Ilower); % lower 2.5% quantile
tauupper(1)=xx(Iupper); % upper 2.5% quantile
% −−−−−−−

ind = randsample(N,N,true,w); % selection
part = part(ind);
for k = 1:n % main loop
    R = A + (B-A).*rand(N,1);
    part = R.*part.*(1-part); % mutation
    w = p(part,Y(k+1)); % weighting
    tau(k + 1) = sum(part.*w)/sum(w); % estimation

    % conf int
    [xx,I]=sort(part); % sort data
    cw=cumsum(w(I))/sum(w); % cumulative normalized weightsum
    % for sorted data
    Ilower=find(cw>=0.025,1); % index for lower 2.5% quantile
    Iupper=find(cw>=0.975,1); % index upper 2.5% quantile
    taulower(k+1)=xx(Ilower); % lower 2.5% quantile
    tauupper(k+1)=xx(Iupper); % upper 2.5% quantile
    % −−−−−−−
    ind = randsample(N,N,true,w); % selection
    part = part(ind);

end

figure(1)
plot(Y, X, '*')
title('Measured data')
xlabel('Y / Observed')
ylabel('X / Hidden')
figure(2)
plot(Y, tau', '*')
title('Estimated data')
xlabel('Y / Observed')
ylabel('\tau_k / Estimated')
figure(3) 
plot(tau, 'b--')
hold on
plot(X', 'ro')
title('Hidden vs Estimated')
legend('\tau_k / Estimated', 'X / Hidden')

figure(4)
plot(X', 'r*')
hold on
plot(taulower, 'b')
hold on
plot(tauupper, 'b')
title('Confidence interval')
legend('X / Hidden', ' Confidence Bounds', 'Location', 'southeast')
hold off

% figure(3)
% plot(X, tau', '*')



