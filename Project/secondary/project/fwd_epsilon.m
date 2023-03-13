function epsilon = fwd_epsilon(data,t,dsigma)
    E = data(1);
    K = data(2);
    n = data(3);
    epsilon = (1/(n+1))*((dsigma/K)^n)*(t.^(n+1)) + (dsigma/E)*t;
end