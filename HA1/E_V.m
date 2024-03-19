function ev = E_V(lambda, k, m)
%Calculates expected value of Weibull distribution
    ev = gamma(1 + m./k).*lambda.^m;
end