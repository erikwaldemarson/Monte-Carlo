function ln_f = ln_f(t, lambda, tau) %f(t | theta, lambda, tau)
    d = length(lambda);
    n_disasters = zeros(1,d)';

    for i = 1:d
        n_disasters(i) = sum(tau >= t(i) & tau < t(i+1)); %number of distasters in interval
    end

    ln_f = 0;
    for i = 1:d
        ln_f = ln_f + log(t(i+1)-t(i)) + n_disasters(i)*log(lambda(i)) - lambda(i)*(t(i+1)-t(i));
    end

   
end
