function epsilon = fwd_epsilon(t,dsigma,E,K,n)
    epsilon = (1/(n+1))*((dsigma/K)^n)*(t.^(n+1)) + (dsigma/E)*t;
end