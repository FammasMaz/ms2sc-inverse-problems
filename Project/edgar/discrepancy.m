function discr = discrepancy (t,dsigma,t_range,epsilon_data,parameters)
    E = parameters(1);
    K = parameters(2);
    n = parameters(3);

    eps = fwd_epsilon(t,dsigma,E,K,n);
    ind=dsearchn(t_range',t');
    discr = 0.5*(eps - epsilon_data(ind)).^2;
end