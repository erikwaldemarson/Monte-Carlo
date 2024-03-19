function ln_f2 = ln_f2(tau, lambda, t) %f(tau | lambda, t)
    d = length(lambda);
    n_disasters = zeros(1,d)';

    for i = 1:d
        n_disasters(i) = sum(tau >= t(i) & tau < t(i+1)); %number of distasters in interval
    end

    ln_f2 = 0;
    for i = 1:d
        ln_f2 = ln_f2 + n_disasters(i)*log(lambda(i)) - lambda(i)*(t(i+1)-t(i));
    end

end